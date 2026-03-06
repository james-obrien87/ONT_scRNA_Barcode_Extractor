#!/usr/bin/env python3
"""
=============================================================================
ONT Single-Cell RNA-seq Barcode Extraction Pipeline
=============================================================================
Author: James O'Brien
        University of Cambridge
        
Description:
    Processes ONT cDNA reads with the following structure (3' end):
        [cDNA insert]-[BC1:8bp]-[Linker1:12bp]-[BC2:8bp]-[Linker2:12bp]-[BC3:8bp]-[UMI:10bp]

    Pipeline stages:
        Stage 0: ONT adapter trimming
        Stage 1: Orientation detection — forward pass, then RC if forward fails
        Stage 2: K-mer anchoring of Linker1 and Linker2
        Stage 3: Levenshtein correction of linker positions
        Stage 4: Anchored barcode + UMI extraction (indel-tolerant window alignment)
        Stage 5: Whitelist correction of BC1, BC2, BC3
        Stage 6: Chimeric read detection (adapter-in-insert + second BC structure)
        Stage 7: cDNA insert extraction + FLAMES-style FASTQ output
        Stage 8: TSV QC summary report

    Read structure (3' end, forward strand):
        [cDNA insert] - [polyA] - [BC1:8bp] - [Linker1:12bp] - [BC2:8bp]
                      - [Linker2:12bp] - [BC3:8bp] - [UMI:10bp]

    Header format (FLAMES-style):
        @<original_read_id>_#<part>_<+/->of<total> CB:Z:... CR:Z:... UR:Z:...
        where + = barcode found on forward strand, - = reverse complement

Usage:
    python extractort.py --config config.yaml --input reads.fastq.gz --output_dir results/

Dependencies:
    pip install pyyaml tqdm
    (edlib is optional — the pipeline uses a built-in Levenshtein implementation
     if edlib is not available, or edlib if it is installed for extra speed)
=============================================================================
"""

import os
import sys
import gzip
import argparse
import logging
import yaml
import multiprocessing
from itertools import islice
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict, List

try:
    import edlib as _edlib
    _HAVE_EDLIB = True
except ImportError:
    _HAVE_EDLIB = False

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------
# Translation table for reverse complementing DNA.
# Handles both cases and preserves N/n (ambiguous bases from basecallers).
_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Translates each base to its complement (A↔T, C↔G) then reverses the
    string. Both upper and lower case are handled; N/n ambiguous bases are
    preserved as N/n. Quality strings are NOT handled here — the caller is
    responsible for reversing the quality string independently (qual[::-1]).

    Args:
        seq: DNA sequence string (may be upper or lower case)

    Returns:
        Reverse-complemented sequence in the same case as input.
    """
    return seq.translate(_RC_TABLE)[::-1]

# ---------------------------------------------------------------------------
# Edit distance and alignment — pure Python fallback + edlib wrapper
#
# Two alignment modes are used throughout the pipeline:
#
#   Global (NW) — used for whitelist barcode correction.
#       Both sequences are aligned end-to-end. Appropriate when comparing
#       a raw barcode against a whitelist entry of the same length.
#
#   Infix / semi-global (HW) — used for adapter trimming, linker finding,
#       and barcode window alignment. The query (adapter/linker/barcode) is
#       found as a substring of the target (read) with free gaps at the
#       target ends. This allows the query to appear anywhere in the target
#       without penalising the unaligned flanking sequence.
#
# When edlib is installed (recommended), both modes delegate to edlib's C
# implementation of the Myers bit-vector algorithm, which is O(NM/w) and
# substantially faster than the pure-Python fallback below.
# ---------------------------------------------------------------------------

def _levenshtein(s: str, t: str, k: int = None) -> int:
    """
    Pure-Python global Levenshtein edit distance between strings s and t.

    Uses the standard two-row rolling DP array (O(min(|s|,|t|)) space).
    Includes an early-exit optimisation: if the row minimum exceeds k at
    any point, the distance is already > k so we return k+1 immediately.
    This makes the function efficient for tight-tolerance whitelist scans
    even without edlib.

    Args:
        s, t: Input strings to compare.
        k:    If provided, exit early and return k+1 when distance > k.

    Returns:
        Integer edit distance, or k+1 if distance exceeds k.
    """
    if s == t:
        return 0
    len_s, len_t = len(s), len(t)
    # Ensure s is the longer string so we allocate the smaller array
    if len_s < len_t:
        s, t, len_s, len_t = t, s, len_t, len_s
    if len_t == 0:
        return len_s

    prev = list(range(len_t + 1))
    for i, cs in enumerate(s):
        curr = [i + 1] + [0] * len_t
        row_min = curr[0]
        for j, ct in enumerate(t):
            curr[j + 1] = min(
                prev[j + 1] + 1,                        # deletion from s
                curr[j]    + 1,                         # insertion into s
                prev[j]    + (0 if cs == ct else 1),    # substitution
            )
            row_min = min(row_min, curr[j + 1])
        # Early exit: if the best possible remaining score already exceeds k,
        # no alignment of the remaining characters can bring it back under k
        if k is not None and row_min > k:
            return k + 1
        prev = curr
    return prev[len_t]


def _align_infix(query: str, target: str, k: int) -> Tuple[int, int, int]:
    """
    Pure-Python infix (semi-global / HW mode) alignment.

    Finds `query` as a substring of `target`, minimising edit distance.
    Gaps at the START and END of the target are free (i.e. the query can
    begin anywhere in the target without penalty). This is equivalent to
    edlib's HW mode and is correct for:
        - Adapter trimming (adapter appears at a variable offset from read end)
        - Linker finding    (linker appears at a roughly known offset in read)
        - Barcode alignment (barcode inside a bounded window)

    Implementation:
        Standard semi-global DP with free first-column initialisation.
        prev[j] represents the cost of aligning query[0..j-1] optimally
        ending at target position i. The last row value gives the minimum
        cost of aligning the full query ending at each target position.

    Args:
        query:  Short sequence to locate (adapter / linker / barcode).
        target: Longer sequence to search in (read or window slice).
        k:      Maximum edit distance; hits with dist > k are ignored.

    Returns:
        (start, end, edit_distance) where end is exclusive.
        Returns (-1, -1, -1) if no hit within k edits is found.
    """
    q_len = len(query)
    t_len = len(target)
    if q_len == 0:
        return 0, 0, 0

    INF = 10**9
    # Initialise: aligning i query chars costs i deletions regardless of target start
    prev = list(range(q_len + 1))

    best_dist = INF
    best_end  = -1

    for j in range(t_len):
        curr = [0] * (q_len + 1)
        curr[0] = 0   # Free gap at target start — no penalty for beginning anywhere
        for i in range(1, q_len + 1):
            curr[i] = min(
                prev[i]     + 1,   # query deletion (advance target, stay in query)
                curr[i - 1] + 1,   # query insertion (advance query, stay in target)
                prev[i - 1] + (0 if query[i - 1] == target[j] else 1),  # match/sub
            )
        # curr[q_len] = cost of aligning the full query ending at target[j]
        if curr[q_len] < best_dist:
            best_dist = curr[q_len]
            best_end  = j + 1   # exclusive end position in target
        prev = curr

    if best_dist > k:
        return -1, -1, -1

    # Recover approximate start position.
    # A full traceback would be O(NM) memory; instead we use an estimate:
    # the alignment spans roughly q_len +/- best_dist characters in the target.
    estimated_start = max(0, best_end - q_len - best_dist)
    return estimated_start, best_end, best_dist


def edit_align_infix(query: str, target: str, k: int) -> Tuple[int, int, int]:
    """
    Public wrapper: find `query` inside `target` (infix/HW mode) within k edits.

    Delegates to edlib if available (fast C implementation), otherwise falls
    back to the pure-Python `_align_infix` above.

    Args:
        query:  Sequence to search for (adapter, linker, or barcode).
        target: Sequence to search within.
        k:      Maximum edit distance threshold.

    Returns:
        (start, end, edit_distance) — end is exclusive.
        Returns (-1, -1, -1) if not found within k edits.
    """
    if _HAVE_EDLIB:
        result = _edlib.align(query, target, mode="HW", task="locations", k=k)
        if result["editDistance"] == -1:
            return -1, -1, -1
        loc = result["locations"][0]
        return loc[0], loc[1] + 1, result["editDistance"]
    else:
        return _align_infix(query, target, k)


def edit_distance(s: str, t: str, k: int = None) -> int:
    """
    Public wrapper: global edit distance between strings s and t.

    Uses edlib (NW mode) if available, otherwise the pure-Python fallback.
    Used for whitelist barcode correction where both strings are the same
    nominal length and we want end-to-end alignment.

    Args:
        s, t: Strings to compare.
        k:    Optional early-exit threshold (returns k+1 if dist > k).

    Returns:
        Integer edit distance, or k+1 / a large number if distance exceeds k.
    """
    if _HAVE_EDLIB:
        result = _edlib.align(s, t, mode="NW", task="distance",
                               k=k if k is not None else -1)
        if result["editDistance"] == -1:
            return (k + 1) if k is not None else 10**9
        return result["editDistance"]
    else:
        return _levenshtein(s, t, k=k)


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    """
    All runtime configuration parameters, parsed from config.yaml.

    Each field maps directly to a YAML key (see load_config for the mapping).
    Default values are safe starting points; most users will need to set at
    minimum the adapter sequences, linker sequences, barcode lengths, and
    whitelist paths.
    """

    # --- Adapter trimming ---
    # ONT adapters are present at both ends of reads in either orientation.
    # trim_adapter() tries the forward and RC of each adapter at each end,
    # so these should always be the canonical forward-strand sequences.
    adapter_5prime: str = ""          # Canonical 5' adapter (forward strand)
    adapter_3prime: str = ""          # Canonical 3' adapter (forward strand)
    adapter_max_errors: int = 2       # Levenshtein tolerance for adapter matching

    # --- Linker sequences and detection tolerances ---
    # Linkers are the fixed sequences flanking BC1/BC2 and BC2/BC3. Their
    # positions serve as the hard structural anchors for barcode extraction.
    linker1_seq: str = ""             # Full Linker1 sequence (e.g. 12 bp)
    linker1_kmer_size: int = 8        # Seed k-mer length for Stage 2 seeding
    linker1_max_levenshtein: int = 2  # Max edit distance for Linker1 acceptance

    linker2_seq: str = ""             # Full Linker2 sequence (e.g. 12 bp)
    linker2_kmer_size: int = 8        # Seed k-mer length for Stage 2 seeding
    linker2_max_levenshtein: int = 2  # Max edit distance for Linker2 acceptance

    # --- Barcode and UMI lengths ---
    # These must match the actual library design. Changing them shifts all
    # window boundaries and anchor offsets throughout the pipeline.
    bc1_length: int = 8               # BC1 length in bp
    bc2_length: int = 8               # BC2 length in bp
    bc3_length: int = 8               # BC3 length in bp
    umi_length: int = 10              # UMI length in bp (anchored to read end)

    # --- Barcode whitelist correction tolerances ---
    # 1 is the recommended starting point. Increasing to 2 risks ambiguous
    # corrections unless the whitelist has high minimum Hamming distance.
    bc1_max_levenshtein: int = 1      # Max edits allowed when correcting BC1
    bc2_max_levenshtein: int = 1      # Max edits allowed when correcting BC2
    bc3_max_levenshtein: int = 1      # Max edits allowed when correcting BC3

    # --- Whitelist file paths ---
    # Plain-text files, one barcode per line. Set to None to skip correction
    # for that barcode (raw extracted sequence is passed through as-is).
    whitelist_bc1: Optional[str] = None
    whitelist_bc2: Optional[str] = None
    whitelist_bc3: Optional[str] = None

    # --- Insert length filter ---
    # Reads whose cDNA insert (everything upstream of the polyA/BC1) is
    # shorter than this after barcode removal are discarded.
    min_insert_length: int = 100      # Minimum cDNA insert length in bp

    # --- PolyA tail detection ---
    # The polyA tail immediately precedes BC1 and acts as the LEFT hard
    # anchor for BC1 extraction. Scanning from the right (toward BC1) finds
    # the boundary closest to BC1 and avoids internal A-runs in the cDNA.
    polya_min_run: int = 6            # Min consecutive A's to call a polyA tail
    polya_max_mismatches: int = 1     # Non-A bases tolerated within the run

    # --- Strand orientation ---
    # ONT reads arrive in either orientation. With try_reverse_complement=True,
    # the pipeline tries the forward strand first, then takes the RC and retries
    # if all barcode stages fail. The successful strand is recorded as '+'/'-'.
    # RC adapter sequences are also tried during adapter trimming (Stage 0)
    # so trimming is always orientation-agnostic regardless of this setting.
    try_reverse_complement: bool = True

    # --- Chimeric read detection ---
    # A read is flagged chimeric if EITHER test passes:
    #   Test A (adapter_check):    an internal adapter is found inside the insert
    #   Test B (structure_check):  a second complete Linker1+Linker2 structure
    #                              is found upstream of the primary barcodes
    # Chimeric reads are written to chimeric_reads.fastq for manual inspection.
    chimeric_adapter_check: bool = True    # Enable Test A
    chimeric_structure_check: bool = True  # Enable Test B
    chimeric_adapter_seq: str = ""         # Adapter to search for internally;
                                           # defaults to adapter_5prime if empty
    chimeric_adapter_max_errors: int = 3   # Levenshtein tolerance for Test A

    # --- Output options ---
    write_failed: bool = True  # Write failed reads to unassigned_reads.fastq


@dataclass
class ReadResult:
    """Result object for a single processed read."""
    read_id: str = ""
    original_seq: str = ""
    original_qual: str = ""

    # After adapter trimming
    trimmed_seq: str = ""
    trimmed_qual: str = ""
    adapter_trimmed: bool = False

    # Strand orientation
    # "+" = barcode found on forward strand of the input read
    # "-" = barcode found only after reverse complementing the read
    strand: str = "+"

    # Linker positions (start index, end index) in the oriented working sequence
    linker1_start: int = -1
    linker1_end: int = -1
    linker1_edit_dist: int = -1

    linker2_start: int = -1
    linker2_end: int = -1
    linker2_edit_dist: int = -1

    # PolyA anchor (end position of polyA tail in working sequence)
    polya_end: int = -1           # -1 = not found; extraction proceeds without polyA anchor

    # Barcode window lengths (for QC — reveals if indels caused non-nominal windows)
    bc1_window_len: int = -1      # actual bp between polyA end and Linker1 start
    bc2_window_len: int = -1      # actual bp between Linker1 end and Linker2 start
    bc3_window_len: int = -1      # actual bp between Linker2 end and UMI region start

    # Extracted barcode strings (raw best-alignment hit from window)
    # On a "-" strand read, these are the RC'd sequences (i.e. forward-strand equivalent)
    bc1_raw: str = ""
    bc2_raw: str = ""
    bc3_raw: str = ""
    umi_raw: str = ""

    # Corrected barcodes (after whitelist lookup; always forward-strand)
    bc1_corrected: str = ""
    bc2_corrected: str = ""
    bc3_corrected: str = ""
    bc1_edit_dist: int = -1
    bc2_edit_dist: int = -1
    bc3_edit_dist: int = -1

    # cDNA insert — always written 5'→3' relative to original mRNA
    # (RC'd if the barcode was found on the "-" strand)
    insert_seq: str = ""
    insert_qual: str = ""

    # Chimeric read fields
    # is_chimeric = True if an internal adapter or second BC structure was found
    is_chimeric: bool = False
    chimeric_reason: str = ""     # "INTERNAL_ADAPTER" | "SECOND_BC_STRUCTURE" | "BOTH"

    # For chimeric sub-reads written from a split: part index and total
    chimeric_part: int = 1        # 1-based index of this sub-read
    chimeric_total: int = 1       # total number of sub-reads from this original read

    # Pass/fail tracking
    passed: bool = False
    fail_reason: str = ""


# ---------------------------------------------------------------------------
# Config loader
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> PipelineConfig:
    """
    Parse a YAML config file and return a populated PipelineConfig dataclass.

    All keys are optional in the YAML — missing keys fall back to the
    PipelineConfig defaults. This makes it safe to add new config sections
    in future without breaking existing config files.

    Args:
        config_path: Path to config.yaml

    Returns:
        PipelineConfig with all fields set from the YAML or defaults.
    """
    with open(config_path) as f:
        raw = yaml.safe_load(f)

    cfg = PipelineConfig()

    # Adapter trimming — forward-strand sequences; RC variants are computed
    # inside trim_adapter() so these should always be the canonical sequences
    a = raw.get("adapter", {})
    cfg.adapter_5prime       = a.get("sequence_5prime", "")
    cfg.adapter_3prime       = a.get("sequence_3prime", "")
    cfg.adapter_max_errors   = int(a.get("max_errors", 2))

    # Linker 1 (between BC1 and BC2)
    l1 = raw.get("linker1", {})
    cfg.linker1_seq             = l1.get("sequence", "").upper()
    cfg.linker1_kmer_size       = int(l1.get("kmer_size", 8))
    cfg.linker1_max_levenshtein = int(l1.get("max_levenshtein", 2))

    # Linker 2 (between BC2 and BC3)
    l2 = raw.get("linker2", {})
    cfg.linker2_seq             = l2.get("sequence", "").upper()
    cfg.linker2_kmer_size       = int(l2.get("kmer_size", 8))
    cfg.linker2_max_levenshtein = int(l2.get("max_levenshtein", 2))

    # Barcode and UMI dimensions + correction tolerances
    bc = raw.get("barcodes", {})
    cfg.bc1_length           = int(bc.get("bc1_length", 8))
    cfg.bc2_length           = int(bc.get("bc2_length", 8))
    cfg.bc3_length           = int(bc.get("bc3_length", 8))
    cfg.umi_length           = int(bc.get("umi_length", 10))
    cfg.bc1_max_levenshtein  = int(bc.get("bc1_max_levenshtein", 1))
    cfg.bc2_max_levenshtein  = int(bc.get("bc2_max_levenshtein", 1))
    cfg.bc3_max_levenshtein  = int(bc.get("bc3_max_levenshtein", 1))
    cfg.whitelist_bc1        = bc.get("whitelist_bc1", None)
    cfg.whitelist_bc2        = bc.get("whitelist_bc2", None)
    cfg.whitelist_bc3        = bc.get("whitelist_bc3", None)

    # Minimum cDNA insert length after barcode removal
    ins = raw.get("insert", {})
    cfg.min_insert_length = int(ins.get("min_length", 100))

    # PolyA tail detection parameters
    pa = raw.get("polya", {})
    cfg.polya_min_run        = int(pa.get("min_run", 6))
    cfg.polya_max_mismatches = int(pa.get("max_mismatches", 1))

    # Strand orientation — whether to try the RC if forward extraction fails
    ori = raw.get("orientation", {})
    cfg.try_reverse_complement = bool(ori.get("try_reverse_complement", True))

    # Chimeric read detection settings
    ch = raw.get("chimeric", {})
    cfg.chimeric_adapter_check    = bool(ch.get("adapter_check", True))
    cfg.chimeric_structure_check  = bool(ch.get("structure_check", True))
    cfg.chimeric_adapter_seq      = ch.get("adapter_seq", "").upper()
    cfg.chimeric_adapter_max_errors = int(ch.get("adapter_max_errors", 3))

    # Output options
    out = raw.get("output", {})
    cfg.write_failed = bool(out.get("write_failed", True))

    log.info("Configuration loaded from %s", config_path)
    return cfg


# ---------------------------------------------------------------------------
# Whitelist loading + barcode correction
# ---------------------------------------------------------------------------

def load_whitelist(path: Optional[str]) -> Optional[set]:
    """
    Load a barcode whitelist from a plain-text file into a Python set.

    The set enables O(1) exact-match lookups, which are done first in
    correct_barcode() before falling through to the slower Levenshtein scan.
    Entries are uppercased on load so comparison is case-insensitive.

    Args:
        path: Path to whitelist file (one barcode per line), or None.

    Returns:
        Set of barcode strings, or None if path is None (no whitelist).
    """
    if path is None:
        return None
    if not os.path.exists(path):
        raise FileNotFoundError(f"Whitelist not found: {path}")
    with open(path) as f:
        wl = set(line.strip().upper() for line in f if line.strip())
    log.info("Loaded whitelist: %s (%d entries)", path, len(wl))
    return wl


def correct_barcode(
    raw_bc: str,
    whitelist: set,
    max_dist: int,
) -> Tuple[str, int]:
    """
    Correct a raw extracted barcode against a whitelist using Levenshtein distance.

    Correction strategy (in order):
        1. Exact match against whitelist — fastest path, distance 0.
        2. Scan all whitelist entries for any within max_dist edits.
        3. If exactly one closest entry → return it with its edit distance.
        4. If two or more entries are tied at the same minimum distance
           (ambiguous) → fail; return ("", distance).
        5. If no entry within max_dist → fail; return ("", -1).

    When no whitelist is provided, the raw barcode is returned unchanged
    with distance 0. This allows the pipeline to run without whitelists
    for exploratory analysis.

    Args:
        raw_bc:    Raw barcode sequence extracted from the read.
        whitelist: Set of known valid barcode sequences (or None).
        max_dist:  Maximum Levenshtein distance for a valid correction.

    Returns:
        (corrected_barcode, edit_distance)
        corrected_barcode is "" if correction failed or was ambiguous.
    """
    if not whitelist:
        return raw_bc, 0

    # Fast exact-match check avoids scanning the whole whitelist for perfect reads
    if raw_bc in whitelist:
        return raw_bc, 0

    best_bc    = ""
    best_dist  = max_dist + 1
    ambiguous  = False

    for wl_bc in whitelist:
        dist = edit_distance(raw_bc, wl_bc, k=max_dist)
        if dist > max_dist:
            continue
        if dist < best_dist:
            # New best candidate — clear any previous ambiguity
            best_dist = dist
            best_bc   = wl_bc
            ambiguous = False
        elif dist == best_dist:
            # Tie with existing best — flag as ambiguous
            if wl_bc != best_bc:
                ambiguous = True

    if ambiguous or best_bc == "":
        return "", best_dist if best_dist <= max_dist else -1
    return best_bc, best_dist


# ---------------------------------------------------------------------------
# Stage 0: Adapter trimming (orientation-aware)
# ---------------------------------------------------------------------------

def trim_adapter(
    seq: str,
    qual: str,
    adapter_5p: str,
    adapter_3p: str,
    max_errors: int,
) -> Tuple[str, str, bool]:
    """
    Trim ONT adapter sequences from both ends of a read, handling both
    forward and reverse-strand reads transparently.

    ONT reads can be sequenced in either orientation.  The adapter that
    appears at the 5' end of a forward-strand read appears as its reverse
    complement at the 3' end of a reverse-strand read (and vice versa).
    To handle this without knowing the read orientation in advance, each
    end is searched with BOTH the canonical adapter and its RC:

        5' end candidates:  adapter_5p  |  RC(adapter_3p)
        3' end candidates:  adapter_3p  |  RC(adapter_5p)

    For each end the best-scoring hit (lowest edit distance) is used.
    This means adapter trimming is fully orientation-agnostic and works
    correctly regardless of which strand the pore sequenced — the Stage 1
    orientation detection then works on a cleanly trimmed sequence.

    Uses semi-global / infix alignment (HW mode): the adapter is found
    as a substring of the read with free gaps at the read ends, which
    is the correct mode for adapter trimming at variable positions.

    Args:
        seq:         Raw read sequence (uppercase)
        qual:        Corresponding quality string
        adapter_5p:  Known 5' adapter sequence (forward strand)
        adapter_3p:  Known 3' adapter sequence (forward strand)
        max_errors:  Maximum Levenshtein distance for a valid adapter hit

    Returns:
        (trimmed_seq, trimmed_qual, was_trimmed)
        was_trimmed is True if at least one adapter end was removed.
    """
    trimmed = False
    start   = 0
    end     = len(seq)

    # Pre-compute reverse complements once
    rc_5p = reverse_complement(adapter_5p) if adapter_5p else ""
    rc_3p = reverse_complement(adapter_3p) if adapter_3p else ""

    # ------------------------------------------------------------------
    # 5' end: try adapter_5p (forward read) and RC(adapter_3p) (RC read)
    # Accept whichever gives the lower edit distance, breaking ties by
    # choosing the hit that removes the most sequence (further trim point).
    # ------------------------------------------------------------------
    best_5p_cut = 0   # exclusive end of the best 5' adapter hit
    best_5p_dist = max_errors + 1

    for candidate in filter(None, [adapter_5p, rc_3p]):
        loc_s, loc_e, dist = edit_align_infix(candidate, seq, max_errors)
        if loc_s == -1:
            continue
        # Prefer lower edit distance; on ties prefer the cut that removes more
        if dist < best_5p_dist or (dist == best_5p_dist and loc_e > best_5p_cut):
            best_5p_dist = dist
            best_5p_cut  = loc_e

    if best_5p_cut > 0:
        start   = best_5p_cut
        trimmed = True

    # ------------------------------------------------------------------
    # 3' end: try adapter_3p (forward read) and RC(adapter_5p) (RC read)
    # Search only within the sequence that remains after 5' trimming.
    # Accept whichever gives the lower edit distance, breaking ties by
    # choosing the hit that removes the most sequence (earlier cut point).
    # ------------------------------------------------------------------
    best_3p_cut  = end   # start of the best 3' adapter hit (exclusive)
    best_3p_dist = max_errors + 1

    if end > start:
        sub = seq[start:end]
        for candidate in filter(None, [adapter_3p, rc_5p]):
            sub_s, sub_e, dist = edit_align_infix(candidate, sub, max_errors)
            if sub_s == -1:
                continue
            abs_cut = start + sub_s   # position in full seq
            if dist < best_3p_dist or (dist == best_3p_dist and abs_cut < best_3p_cut):
                best_3p_dist = dist
                best_3p_cut  = abs_cut

    if best_3p_cut < end:
        end     = best_3p_cut
        trimmed = True

    return seq[start:end], qual[start:end], trimmed


# ---------------------------------------------------------------------------
# Stage 2: K-mer anchoring + Levenshtein linker search
#
# Linker detection is a two-pass process designed to be both fast and
# tolerant of ONT basecalling errors:
#
#   Pass 1 — K-mer seeding (O(n)):
#       Build a k-mer index of the known linker sequence. Slide over the
#       restricted 3' search region of the read. Wherever a read k-mer
#       matches a linker k-mer, record the implied linker start position
#       as a candidate. This avoids full alignment of the entire read.
#
#   Pass 2 — Levenshtein verification (O(candidates × linker_len)):
#       For each candidate start, extract a window of width
#       (linker_len ± max_edits) and compute the infix edit distance
#       between the known linker and the window. Accept the best hit
#       within the configured tolerance.
#
# The search region for Linker1 is restricted to the expected 3' barcode
# region (last bc1+L1+bc2+L2+bc3+UMI + buffer bases of the read) to
# avoid spurious matches in A/T-rich regions of the cDNA body.
# Linker2 is searched only after Linker1's end position, enforcing the
# correct structural order.
# ---------------------------------------------------------------------------

def build_kmer_index(sequence: str, k: int) -> Dict[str, List[int]]:
    """
    Build a k-mer → [positions] index for a given sequence.

    For a linker of length L with k-mer size k, this creates L-k+1 entries.
    Each entry maps a k-mer substring to the list of positions within the
    linker where that k-mer occurs (usually just one position per k-mer for
    a non-repetitive linker sequence).

    Args:
        sequence: The sequence to index (linker sequence).
        k:        K-mer length.

    Returns:
        Dict mapping each k-mer to a list of its start positions in sequence.
    """
    index = defaultdict(list)
    for i in range(len(sequence) - k + 1):
        index[sequence[i:i+k]].append(i)
    return index


def find_linker(
    seq: str,
    linker: str,
    kmer_size: int,
    max_levenshtein: int,
    search_start: int = 0,
    search_end: Optional[int] = None,
) -> Tuple[int, int, int]:
    """
    Locate a linker sequence within a read using k-mer seeding + Levenshtein
    verification (see section header above for the full algorithm description).

    The search is restricted to seq[search_start:search_end] to avoid
    spurious matches in the cDNA body. Returned positions are always in
    full-sequence coordinates (not relative to search_start).

    Args:
        seq:             Full read sequence (after adapter trimming).
        linker:          Known linker sequence to find.
        kmer_size:       K-mer seed length (smaller = more sensitive, more noise).
        max_levenshtein: Maximum edit distance for an accepted hit.
        search_start:    Start of the restricted search window in seq.
        search_end:      End of the restricted search window in seq (exclusive).
                         Defaults to len(seq) if None.

    Returns:
        (linker_start, linker_end, edit_distance) in full seq coordinates.
        linker_end is exclusive (Python slice convention).
        Returns (-1, -1, -1) if no hit within tolerance.
    """
    if search_end is None:
        search_end = len(seq)

    sub_seq    = seq[search_start:search_end]
    linker_len = len(linker)

    # Build k-mer index of the linker sequence
    linker_kmers = build_kmer_index(linker, kmer_size)

    # Slide over the search region and collect candidate linker start positions
    # by finding exact k-mer matches. Each hit implies an estimated linker start.
    candidate_positions = set()
    for i in range(len(sub_seq) - kmer_size + 1):
        kmer = sub_seq[i:i+kmer_size]
        if kmer in linker_kmers:
            for linker_pos in linker_kmers[kmer]:
                # If read[i] matches linker[linker_pos], linker starts at i - linker_pos
                est_start = i - linker_pos
                candidate_positions.add(est_start)

    if not candidate_positions:
        return -1, -1, -1

    # Verify each candidate with a full infix alignment in a padded window.
    # The window is padded by ±max_levenshtein to accommodate indels that
    # shift the linker by up to max_levenshtein bases relative to the seed estimate.
    best_start = -1
    best_end   = -1
    best_dist  = max_levenshtein + 1
    slack      = max_levenshtein

    for cand in candidate_positions:
        win_start = max(0, cand - slack)
        win_end   = min(len(sub_seq), cand + linker_len + slack)
        window    = sub_seq[win_start:win_end]

        loc_start, loc_end, dist = edit_align_infix(linker, window, max_levenshtein)
        if loc_start == -1:
            continue
        if dist < best_dist:
            best_dist  = dist
            # Convert window-relative positions back to full seq coordinates
            best_start = search_start + win_start + loc_start
            best_end   = search_start + win_start + loc_end  # already exclusive

    if best_start == -1:
        return -1, -1, -1
    return best_start, best_end, best_dist


# ---------------------------------------------------------------------------
# Stage 3: PolyA tail detection (LEFT anchor for BC1)
#
# The polyA tail is the last structural feature of the cDNA molecule and
# immediately precedes BC1. Confirming its position serves two purposes:
#
#   1. LEFT ANCHOR for BC1 extraction:
#      The polyA end is the hard left boundary of the BC1 window. Combined
#      with the Linker1 start position on the right, BC1 is fully enclosed
#      between two confirmed anchors, making the window robust to ±1–2 bp
#      indels in BC1 without any fixed-width slicing.
#
#   2. Insert boundary:
#      The polyA run is trimmed away when writing the cDNA insert, giving
#      a clean poly-A-free insert sequence for downstream alignment.
#
# Search direction — right-to-left (scan from the 3' side of the insert):
# The polyA tail abuts BC1. Scanning from the right finds the run that is
# closest to BC1 and avoids matching internal A-rich stretches in the cDNA.
# ---------------------------------------------------------------------------

def find_polya_end(
    seq: str,
    search_start: int,
    search_end: int,
    min_run: int,
    max_mismatches: int,
) -> int:
    """
    Locate the end position of a polyA tail in seq[search_start:search_end].

    Scans right-to-left within the search window looking for the rightmost
    run of at least min_run consecutive A's, tolerating up to max_mismatches
    non-A bases within the run. Tolerance is necessary because ONT basecallers
    frequently miscall individual bases within homopolymer runs.

    Once the rightmost qualifying run is found, the function attempts to
    extend it rightward as far as A's continue (or until mismatches are
    exhausted), giving the tightest possible boundary next to BC1.

    Args:
        seq:            Full read sequence (forward strand, already oriented).
        search_start:   Left boundary of the search region.
        search_end:     Right boundary of the search region (exclusive).
                        Should be set to approximately Linker1_start - BC1_length
                        so the search stays upstream of BC1.
        min_run:        Minimum consecutive A's (within mismatch budget) to
                        qualify as a polyA tail. Recommended: 6–8.
        max_mismatches: Maximum non-A bases tolerated within the run.
                        Recommended: 1 for ONT noise tolerance.

    Returns:
        Position of the first base AFTER the polyA run in full seq coordinates.
        This is the expected left boundary of BC1.
        Returns -1 if no qualifying run is found; extraction then falls back
        to a linker-anchored estimate with a small slack window.
    """
    sub     = seq[search_start:search_end]
    sub_len = len(sub)

    best_end = -1  # End of the best polyA run found (exclusive, sub-coords)

    # Slide a window of exactly min_run from right to left.
    # The first qualifying window found is the rightmost one — closest to BC1.
    i = sub_len - min_run
    while i >= 0:
        window     = sub[i:i + min_run]
        mismatches = sum(1 for b in window if b != 'A')
        if mismatches <= max_mismatches:
            # Found a qualifying seed run — try to extend its right boundary
            run_end = i + min_run
            while run_end < sub_len and (
                sub[run_end] == 'A' or
                sum(1 for b in sub[i:run_end + 1] if b != 'A') <= max_mismatches
            ):
                run_end += 1
            best_end = run_end
            break   # Take the rightmost qualifying run; don't search further left
        i -= 1

    if best_end == -1:
        return -1
    # Convert from sub-sequence coordinates to full sequence coordinates
    return search_start + best_end


# ---------------------------------------------------------------------------
# Stage 4: Anchored barcode + UMI extraction
#
# Core design principle: EVERY barcode window has confirmed anchors on BOTH
# sides before any extraction is attempted. This eliminates fixed-width
# slicing and makes the pipeline robust to ±1–2 bp indels in any barcode.
#
# Window boundaries:
#
#   BC1  │ LEFT:  polya_end  (or estimated if polyA not found)
#        │ RIGHT: linker1_start
#        │
#   BC2  │ LEFT:  linker1_end
#        │ RIGHT: linker2_start
#        │
#   BC3  │ LEFT:  linker2_end
#        │ RIGHT: read_end − umi_length
#        │
#   UMI  │ LEFT:  read_end − umi_length   (anchored to read end)
#        │ RIGHT: read_end
#
# Why the UMI is anchored to the read end, not to BC3:
#   The polyA-capture protocol captures reads at a fixed 3' end. The UMI
#   therefore always occupies the last umi_length bases regardless of any
#   indels in BC3. This makes UMI extraction completely independent of BC3.
#
# Within each window, align_barcode_to_window() aligns all whitelist
# candidates into the window in HW (infix) mode. An 8bp barcode with a
# 1bp deletion will span only 7bp of the window — the alignment finds the
# best match regardless of the barcode's actual length in the read.
# ---------------------------------------------------------------------------

def align_barcode_to_window(
    window: str,
    whitelist: Optional[set],
    expected_len: int,
    max_levenshtein: int,
) -> Tuple[str, str, int, int, int]:
    """
    Extract and correct a single barcode from a bounded sequence window.

    The window spans from one confirmed structural anchor to the next, so
    any barcode indel is fully contained within the window. Alignment is
    done in HW (infix) mode: the barcode candidate is found as a substring
    of the window, so flanking anchor bases at the window edge (e.g. the
    last A of the polyA tail) do not penalise the alignment score.

    If a whitelist is provided, every entry is aligned into the window and
    the entry with the lowest edit distance is accepted. Ambiguous ties
    (two different whitelist entries at the same distance) cause a failure
    because a unique cell cannot be assigned.

    If no whitelist is provided, the central expected_len bases of the
    window are extracted as the raw barcode (no correction is applied).
    This is useful for UMI-style extraction or exploratory runs.

    Args:
        window:          Sequence between confirmed left and right anchors.
        whitelist:       Set of known barcode sequences, or None.
        expected_len:    Nominal barcode length in bp (e.g. 8).
        max_levenshtein: Maximum edit distance for acceptance.

    Returns:
        (raw_extracted, corrected, edit_distance, align_start, align_end)
        - raw_extracted: the window substring that best aligned to the winner
        - corrected:     the whitelist entry (or raw if no whitelist); "" on failure
        - edit_distance: Levenshtein distance of the winning alignment
        - align_start/end: position within `window` of the best alignment
        All fields are -1 / "" on failure.
    """
    if not window:
        return "", "", -1, -1, -1

    if whitelist:
        best_bc    = ""
        best_raw   = ""
        best_dist  = max_levenshtein + 1
        best_start = -1
        best_end   = -1
        ambiguous  = False

        for wl_bc in whitelist:
            loc_start, loc_end, dist = edit_align_infix(wl_bc, window, max_levenshtein)
            if loc_start == -1:
                continue
            if dist < best_dist:
                # New best — clear any previous ambiguity flag
                best_dist  = dist
                best_bc    = wl_bc
                best_raw   = window[loc_start:loc_end]
                best_start = loc_start
                best_end   = loc_end
                ambiguous  = False
            elif dist == best_dist and wl_bc != best_bc:
                # Same distance, different barcode — correction is ambiguous
                ambiguous = True

        if ambiguous or not best_bc:
            return "", "", -1, -1, -1
        return best_raw, best_bc, best_dist, best_start, best_end

    else:
        # No whitelist — extract the central expected_len bases of the window.
        # This avoids any anchor bleed (e.g. polyA A's at the left edge)
        # by centering rather than taking from the left.
        mid   = (len(window) - expected_len) // 2
        start = max(0, mid)
        end   = min(len(window), start + expected_len)
        raw   = window[start:end]
        return raw, raw, 0, start, end


def extract_barcodes_anchored(
    seq: str,
    polya_end: int,
    linker1_start: int,
    linker1_end: int,
    linker2_start: int,
    linker2_end: int,
    read_end: int,
    cfg: PipelineConfig,
    wl_bc1: Optional[set],
    wl_bc2: Optional[set],
    wl_bc3: Optional[set],
) -> Tuple[str, str, int, str, str, int, str, str, int, str, int, int, int]:
    """
    Extract BC1, BC2, BC3, and UMI using dual hard anchors at every boundary.

    This is the central extraction function. All four components are extracted
    from windows whose both boundaries are independently confirmed structural
    positions. The table below summarises each window:

        Component  Left anchor           Right anchor
        ─────────  ────────────────────  ───────────────────────────────
        BC1        polya_end             linker1_start
        BC2        linker1_end           linker2_start
        BC3        linker2_end           read_end − umi_length
        UMI        read_end − umi_length read_end   (always fixed)

    BC1 left anchor note:
        If polyA is detected, bc1_left is set to polya_end − bc1_max_levenshtein.
        The small leftward extension ensures a BC1 with a 1bp insertion (which
        shifts the barcode into the polyA region) is still fully captured in
        the window. The extra polyA bases at the window edge don't affect the
        alignment because HW mode allows free gaps in the target flanks.

        If polyA is NOT detected, bc1_left is estimated as
        linker1_start − bc1_length − slack. This is less precise but still
        functional — the window is wider and the alignment remains correct.

    Args:
        seq:            Read sequence in the orientation being processed.
        polya_end:      Position of the first base after the polyA tail,
                        or -1 if polyA was not detected.
        linker1_start/end: Confirmed Linker1 boundaries (end is exclusive).
        linker2_start/end: Confirmed Linker2 boundaries (end is exclusive).
        read_end:       Length of seq (used to anchor BC3 right and UMI).
        cfg:            PipelineConfig (barcode lengths, tolerances).
        wl_bc1/2/3:     Whitelists (or None) for each barcode.

    Returns:
        Flat tuple of 13 values:
            bc1_raw, bc1_corr, bc1_ed,
            bc2_raw, bc2_corr, bc2_ed,
            bc3_raw, bc3_corr, bc3_ed,
            umi_raw,
            bc1_window_len, bc2_window_len, bc3_window_len
        *_raw   = sequence extracted from the window (before correction)
        *_corr  = whitelist-corrected barcode ("" on failure)
        *_ed    = edit distance to whitelist entry (-1 on failure)
        *_window_len = actual window width in bp (useful for QC)
    """
    # --- BC1 window ---
    # Left anchor: polya_end (preferred) or a linker1-relative estimate.
    # We extend left by bc1_max_levenshtein so that a BC1 with a 1bp insertion
    # (which shifts the barcode into the polyA region by 1 base) is still
    # captured within the window.
    indel_slack = cfg.bc1_max_levenshtein
    slack       = cfg.bc1_max_levenshtein + 1   # slightly wider fallback slack
    if polya_end != -1:
        bc1_left = max(0, polya_end - indel_slack)
    else:
        # No polyA — use a wider window anchored from the left of Linker1
        bc1_left = max(0, linker1_start - cfg.bc1_length - slack)

    bc1_window = seq[bc1_left:linker1_start]

    bc1_raw, bc1_corr, bc1_ed, bc1_aln_s, bc1_aln_e = align_barcode_to_window(
        bc1_window, wl_bc1, cfg.bc1_length, cfg.bc1_max_levenshtein
    )

    # --- BC2 window ---
    # Both boundaries are confirmed linker positions — the tightest possible window.
    bc2_window = seq[linker1_end:linker2_start]

    bc2_raw, bc2_corr, bc2_ed, bc2_aln_s, bc2_aln_e = align_barcode_to_window(
        bc2_window, wl_bc2, cfg.bc2_length, cfg.bc2_max_levenshtein
    )

    # --- BC3 window ---
    # Right anchor: read_end − umi_length. The UMI always occupies the last
    # umi_length bases of the read (anchored to the 3' captured end), so BC3
    # ends wherever the UMI begins, regardless of any indel in BC3 itself.
    bc3_right  = read_end - cfg.umi_length
    bc3_window = seq[linker2_end:bc3_right]

    bc3_raw, bc3_corr, bc3_ed, bc3_aln_s, bc3_aln_e = align_barcode_to_window(
        bc3_window, wl_bc3, cfg.bc3_length, cfg.bc3_max_levenshtein
    )

    # --- UMI ---
    # Always the last umi_length bases of the read. This is completely
    # independent of all barcode positions — even a large BC3 indel leaves
    # the UMI position unchanged.
    umi_raw = seq[read_end - cfg.umi_length:read_end]

    bc1_window_len = len(bc1_window)
    bc2_window_len = len(bc2_window)
    bc3_window_len = len(bc3_window)

    return (
        bc1_raw, bc1_corr, bc1_ed,
        bc2_raw, bc2_corr, bc2_ed,
        bc3_raw, bc3_corr, bc3_ed,
        umi_raw,
        bc1_window_len, bc2_window_len, bc3_window_len,
    )


# ---------------------------------------------------------------------------
# Stage 6: Chimeric read detection
# ---------------------------------------------------------------------------

def detect_chimera(
    insert_seq: str,
    cfg: PipelineConfig,
) -> Tuple[bool, str]:
    """
    Inspect the putative cDNA insert for evidence of a chimeric read.

    Two independent tests are run (both enabled by default, each configurable):

    Test A — Internal adapter check:
        Search for the ONT adapter sequence within the insert using infix
        alignment.  An adapter inside the insert strongly indicates that two
        reads were ligated end-to-end during library preparation.
        The adapter sequence used is cfg.chimeric_adapter_seq if set,
        otherwise cfg.adapter_5prime as a fallback.

    Test B — Second barcode structure check:
        Search for Linker1 AND Linker2 within the insert.  If both are found
        in the correct order with sufficient separation to accommodate BC1,
        BC2, BC3, and UMI, a complete second barcode structure exists inside
        the insert, confirming a chimera.

    Args:
        insert_seq: Putative cDNA insert (everything upstream of the polyA/BC1)
        cfg:        PipelineConfig

    Returns:
        (is_chimeric, reason_string)
        reason_string is one of: "", "INTERNAL_ADAPTER", "SECOND_BC_STRUCTURE", "BOTH"
    """
    reasons = []

    # --- Test A: internal adapter ---
    if cfg.chimeric_adapter_check:
        adapter_to_check = cfg.chimeric_adapter_seq or cfg.adapter_5prime
        if adapter_to_check and len(insert_seq) >= len(adapter_to_check):
            loc_s, loc_e, dist = edit_align_infix(
                adapter_to_check,
                insert_seq,
                cfg.chimeric_adapter_max_errors,
            )
            if loc_s != -1:
                reasons.append("INTERNAL_ADAPTER")

    # --- Test B: second barcode structure ---
    if cfg.chimeric_structure_check and len(insert_seq) > 0:
        min_second_structure = (
            cfg.bc1_length + len(cfg.linker1_seq) +
            cfg.bc2_length + len(cfg.linker2_seq) +
            cfg.bc3_length + cfg.umi_length
        )
        if len(insert_seq) >= min_second_structure:
            # Search for Linker1 anywhere in the insert
            l1s, l1e, _ = find_linker(
                insert_seq,
                cfg.linker1_seq,
                cfg.linker1_kmer_size,
                cfg.linker1_max_levenshtein,
            )
            if l1s != -1:
                # Then Linker2 must follow Linker1
                l2s, l2e, _ = find_linker(
                    insert_seq,
                    cfg.linker2_seq,
                    cfg.linker2_kmer_size,
                    cfg.linker2_max_levenshtein,
                    search_start=l1e,
                )
                if l2s != -1:
                    reasons.append("SECOND_BC_STRUCTURE")

    if not reasons:
        return False, ""
    return True, "_AND_".join(reasons)


# ---------------------------------------------------------------------------
# FASTQ record formatting — FLAMES-style header
# ---------------------------------------------------------------------------

def _build_header(
    original_read_id: str,
    strand: str,
    part: int,
    total_parts: int,
    bc1_corr: str,
    bc2_corr: str,
    bc3_corr: str,
    bc1_raw: str,
    bc2_raw: str,
    bc3_raw: str,
    umi_raw: str,
    bc1_ed: int,
    bc2_ed: int,
    bc3_ed: int,
    l1_ed: int,
    l2_ed: int,
    is_chimeric: bool,
    chimeric_reason: str,
) -> str:
    """
    Build a FLAMES-style FASTQ read header.

    Format:
        @<read_id>_#<part>_<strand>of<total> CB:Z:<BC1BC2BC3> CR:Z:<BC1_BC2_BC3>
          UR:Z:<UMI> BC1_ED:i:<d> BC2_ED:i:<d> BC3_ED:i:<d>
          L1_ED:i:<d> L2_ED:i:<d> [XC:Z:CHIMERIC_<reason>]

    Fields:
        _#<part>_<strand>of<total>
            part   = 1-based index of this sub-read within the original read
            strand = "+" forward, "-" reverse complement
            total  = total sub-read count from this original read (1 for non-chimeric)

        CB:Z:  Corrected concatenated cell barcode (BC1+BC2+BC3, forward strand)
        CR:Z:  Raw extracted barcodes (BC1_BC2_BC3, underscore-separated, forward strand)
        UR:Z:  Raw UMI sequence (forward strand)
        BC*_ED  Levenshtein distance to whitelist entry for each barcode
        L*_ED   Levenshtein distance for each linker match
        XC:Z:  Chimeric flag (only present if is_chimeric=True)
    """
    cell_barcode = bc1_corr + bc2_corr + bc3_corr
    raw_barcode  = f"{bc1_raw}_{bc2_raw}_{bc3_raw}"

    name = f"@{original_read_id}_#{part}_{strand}of{total_parts}"

    tags = (
        f"CB:Z:{cell_barcode} "
        f"CR:Z:{raw_barcode} "
        f"UR:Z:{umi_raw} "
        f"BC1_ED:i:{bc1_ed} "
        f"BC2_ED:i:{bc2_ed} "
        f"BC3_ED:i:{bc3_ed} "
        f"L1_ED:i:{l1_ed} "
        f"L2_ED:i:{l2_ed}"
    )
    if is_chimeric:
        tags += f" XC:Z:CHIMERIC_{chimeric_reason}"

    return f"{name} {tags}"


def format_fastq_record(result: ReadResult) -> str:
    """
    Format a passing ReadResult as a FLAMES-style FASTQ record.
    The insert sequence and quality are already in forward-strand orientation.
    """
    header = _build_header(
        original_read_id = result.read_id,
        strand           = result.strand,
        part             = result.chimeric_part,
        total_parts      = result.chimeric_total,
        bc1_corr         = result.bc1_corrected,
        bc2_corr         = result.bc2_corrected,
        bc3_corr         = result.bc3_corrected,
        bc1_raw          = result.bc1_raw,
        bc2_raw          = result.bc2_raw,
        bc3_raw          = result.bc3_raw,
        umi_raw          = result.umi_raw,
        bc1_ed           = result.bc1_edit_dist,
        bc2_ed           = result.bc2_edit_dist,
        bc3_ed           = result.bc3_edit_dist,
        l1_ed            = result.linker1_edit_dist,
        l2_ed            = result.linker2_edit_dist,
        is_chimeric      = result.is_chimeric,
        chimeric_reason  = result.chimeric_reason,
    )
    return f"{header}\n{result.insert_seq}\n+\n{result.insert_qual}\n"


def format_failed_fastq_record(result: ReadResult) -> str:
    """
    Format a failed read for the unassigned FASTQ.
    Uses the trimmed sequence (or original if trimming also failed).
    Embeds fail reason and strand attempt in header.
    """
    seq  = result.trimmed_seq if result.trimmed_seq else result.original_seq
    qual = result.trimmed_qual if result.trimmed_qual else result.original_qual
    header = (
        f"@{result.read_id}_#1_{result.strand}of1 "
        f"FAIL_REASON:{result.fail_reason}"
    )
    return f"{header}\n{seq}\n+\n{qual}\n"


# ---------------------------------------------------------------------------
# Core per-read processing — orientation-aware, chimeric-aware
# ---------------------------------------------------------------------------

def _attempt_extraction(
    working_seq: str,
    working_qual: str,
    strand: str,
    cfg: PipelineConfig,
    wl_bc1: Optional[set],
    wl_bc2: Optional[set],
    wl_bc3: Optional[set],
) -> Optional[dict]:
    """
    Attempt full barcode extraction (Stages 2–5) from a single oriented sequence.

    This function is called by process_read() up to twice: once for the forward
    strand and, if that fails, once for the reverse complement. Separating the
    orientation-independent extraction logic here avoids code duplication and
    keeps process_read() focused on orchestration.

    On success, returns a dict whose keys map directly onto ReadResult fields.
    On any failure (linker not found, barcode correction failed, insert too short),
    returns either None (linker not found — no dict available) or a dict with a
    "fail" key containing the failure reason string.

    Args:
        working_seq:  Read sequence in the orientation being attempted.
        working_qual: Corresponding quality string (already reversed if strand=="-").
        strand:       "+" for forward, "-" for reverse complement.
        cfg:          PipelineConfig.
        wl_bc1/2/3:  Whitelists (or None) for barcode correction.

    Returns:
        dict on success or partial failure (with "fail" key), or None if
        Linker1 was not found (indicating the orientation is likely wrong).
    """
    # ---- Stage 2: Linker 1 search ----------------------------------------
    # Restrict the search to the expected 3' barcode region to avoid spurious
    # matches against A/T-rich regions of the cDNA body. The search window
    # covers: BC1 + L1 + BC2 + L2 + BC3 + UMI + a small buffer for indels.
    barcode_region_len = (
        cfg.bc1_length + len(cfg.linker1_seq) +
        cfg.bc2_length + len(cfg.linker2_seq) +
        cfg.bc3_length + cfg.umi_length
    )
    search_buffer = 20   # Extra bases of slack for read-level indels / soft-clipping
    search_start  = max(0, len(working_seq) - barcode_region_len - search_buffer)

    l1_start, l1_end, l1_dist = find_linker(
        working_seq, cfg.linker1_seq,
        cfg.linker1_kmer_size, cfg.linker1_max_levenshtein,
        search_start=search_start,
    )
    if l1_start == -1:
        # No Linker1 hit — return None to signal a clean orientation mismatch
        # (as opposed to a barcode-level failure which returns a dict with "fail")
        return None

    # ---- Stage 2: Linker 2 search ----------------------------------------
    # Search begins immediately after Linker1 ends to enforce correct order
    # and prevent Linker2 from matching in the BC1 or cDNA regions.
    l2_start, l2_end, l2_dist = find_linker(
        working_seq, cfg.linker2_seq,
        cfg.linker2_kmer_size, cfg.linker2_max_levenshtein,
        search_start=l1_end,
    )
    if l2_start == -1:
        return None

    # ---- Stage 3: PolyA tail detection (left anchor for BC1) -------------
    # Set the search region to end just before the expected BC1 start, so the
    # polyA scanner doesn't accidentally match As within BC1 itself.
    polya_search_end   = l1_start - cfg.bc1_length
    polya_search_start = max(0, polya_search_end - cfg.polya_min_run - 10)

    polya_end = find_polya_end(
        working_seq,
        search_start=polya_search_start,
        search_end=max(polya_search_start, polya_search_end),
        min_run=cfg.polya_min_run,
        max_mismatches=cfg.polya_max_mismatches,
    )
    # polya_end == -1 is valid — extraction degrades gracefully to a
    # linker-anchored estimate; the read is not failed just for missing polyA.

    # ---- Stage 4: Anchored barcode + UMI extraction ----------------------
    read_end = len(working_seq)
    (
        bc1_raw, bc1_corr, bc1_ed,
        bc2_raw, bc2_corr, bc2_ed,
        bc3_raw, bc3_corr, bc3_ed,
        umi_raw,
        bc1_wlen, bc2_wlen, bc3_wlen,
    ) = extract_barcodes_anchored(
        working_seq,
        polya_end,
        l1_start, l1_end,
        l2_start, l2_end,
        read_end,
        cfg,
        wl_bc1, wl_bc2, wl_bc3,
    )

    # ---- Stage 5: Validate barcode correction outcomes -------------------
    # Each barcode must have been corrected to a unique whitelist entry.
    # Returning a dict with "fail" (rather than None) preserves the partial
    # raw barcode information for the unassigned FASTQ header.
    if not bc1_corr:
        return {"fail": "BC1_CORRECTION_FAILED",
                "bc1_raw": bc1_raw, "bc2_raw": bc2_raw, "bc3_raw": bc3_raw}
    if not bc2_corr:
        return {"fail": "BC2_CORRECTION_FAILED",
                "bc1_raw": bc1_raw, "bc2_raw": bc2_raw, "bc3_raw": bc3_raw}
    if not bc3_corr:
        return {"fail": "BC3_CORRECTION_FAILED",
                "bc1_raw": bc1_raw, "bc2_raw": bc2_raw, "bc3_raw": bc3_raw}
    if len(umi_raw) != cfg.umi_length:
        # This is an edge case — should only happen if the read is very short
        return {"fail": f"UMI_EXTRACTION_FAILED_len{len(umi_raw)}"}

    # ---- Stage 5: Extract cDNA insert ------------------------------------
    # The insert is everything upstream of the polyA/BC1 region.
    # If polyA was found, walk left from polya_end to find the true start
    # of the polyA run (excluding the run itself from the insert).
    # If polyA was not found, use BC1's window left edge as the boundary.
    if polya_end != -1:
        pa_start   = polya_end
        mismatches = 0
        while pa_start > 0:
            base = working_seq[pa_start - 1]
            if base != 'A':
                mismatches += 1
                if mismatches > cfg.polya_max_mismatches:
                    break
            pa_start -= 1
        insert_end = pa_start
    else:
        # No polyA — estimate insert end from BC1 window width
        insert_end = max(0, l1_start - bc1_wlen)

    insert_seq  = working_seq[:insert_end]
    insert_qual = working_qual[:insert_end]

    if len(insert_seq) < cfg.min_insert_length:
        return {"fail": f"INSERT_TOO_SHORT_{len(insert_seq)}bp"}

    # All stages passed — return a complete extraction result dict
    return {
        "strand":            strand,
        "linker1_start":     l1_start,
        "linker1_end":       l1_end,
        "linker1_edit_dist": l1_dist,
        "linker2_start":     l2_start,
        "linker2_end":       l2_end,
        "linker2_edit_dist": l2_dist,
        "polya_end":         polya_end,
        "bc1_raw":           bc1_raw,
        "bc2_raw":           bc2_raw,
        "bc3_raw":           bc3_raw,
        "umi_raw":           umi_raw,
        "bc1_corrected":     bc1_corr,
        "bc2_corrected":     bc2_corr,
        "bc3_corrected":     bc3_corr,
        "bc1_edit_dist":     bc1_ed,
        "bc2_edit_dist":     bc2_ed,
        "bc3_edit_dist":     bc3_ed,
        "bc1_window_len":    bc1_wlen,
        "bc2_window_len":    bc2_wlen,
        "bc3_window_len":    bc3_wlen,
        "insert_seq":        insert_seq,
        "insert_qual":       insert_qual,
    }


def process_read(args) -> List[ReadResult]:
    """
    Process a single read through all pipeline stages and return results.

    This is the top-level per-read function called by the multiprocessing pool.
    It orchestrates the full pipeline from adapter trimming through chimeric
    detection and returns a LIST of ReadResult objects so that future support
    for splitting chimeric reads into multiple sub-reads is straightforward.
    For now the list always contains exactly one element.

    Orientation strategy (Stage 1):
        1. Try forward strand — call _attempt_extraction() on the trimmed read.
        2. If forward fails at any stage, reverse complement the trimmed sequence
           and quality string (qual is simply reversed, not RC'd) and retry.
        3. Record the successful strand as "+" or "-" on the ReadResult.
        4. If both fail, record the forward failure reason (primary attempt).

    Chimeric detection (Stage 6):
        After successful barcode extraction, inspect the confirmed insert for
        evidence of a chimeric read. Chimeric reads still pass barcode extraction
        but are routed to chimeric_reads.fastq rather than assigned_reads.fastq.

    Args:
        args: Tuple of (read_id, seq, qual, cfg, wl_bc1, wl_bc2, wl_bc3).
              Packed as a tuple for compatibility with multiprocessing.Pool.map().

    Returns:
        List[ReadResult] — always length 1 in the current implementation.
    """
    read_id, seq, qual, cfg, wl_bc1, wl_bc2, wl_bc3 = args
    r = ReadResult(read_id=read_id, original_seq=seq, original_qual=qual)

    # ------------------------------------------------------------------
    # Stage 0: Adapter trimming
    # Orientation-aware — tries both the canonical adapter and its RC at
    # each read end. See trim_adapter() for the full logic.
    # ------------------------------------------------------------------
    trimmed_seq, trimmed_qual, was_trimmed = trim_adapter(
        seq, qual,
        cfg.adapter_5prime,
        cfg.adapter_3prime,
        cfg.adapter_max_errors,
    )
    r.trimmed_seq     = trimmed_seq
    r.trimmed_qual    = trimmed_qual
    r.adapter_trimmed = was_trimmed

    # Fail fast if the trimmed read is too short to possibly contain the full
    # barcode structure plus the minimum insert length
    min_read_len = (
        cfg.min_insert_length + cfg.bc1_length + len(cfg.linker1_seq) +
        cfg.bc2_length + len(cfg.linker2_seq) + cfg.bc3_length + cfg.umi_length
    )
    if len(trimmed_seq) < min_read_len:
        r.fail_reason = "READ_TOO_SHORT_AFTER_TRIMMING"
        return [r]

    # ------------------------------------------------------------------
    # Stage 1: Orientation detection — forward attempt first
    # _attempt_extraction() returns None if Linker1 was not found at all
    # (clean orientation mismatch), or a dict with "fail" if a later stage
    # failed (barcode correction, insert length, etc.).
    # ------------------------------------------------------------------
    extraction = _attempt_extraction(
        trimmed_seq, trimmed_qual, "+",
        cfg, wl_bc1, wl_bc2, wl_bc3,
    )

    if extraction is None or "fail" in extraction:
        # Forward pass failed — save the deepest failure reason reached
        fwd_fail = (extraction or {}).get("fail", "LINKER1_NOT_FOUND")

        if cfg.try_reverse_complement:
            # Reverse complement the sequence; reverse (not RC) the quality string
            rc_seq  = reverse_complement(trimmed_seq)
            rc_qual = trimmed_qual[::-1]

            extraction = _attempt_extraction(
                rc_seq, rc_qual, "-",
                cfg, wl_bc1, wl_bc2, wl_bc3,
            )
            if extraction is None or "fail" in extraction:
                # Both orientations failed — report the forward failure as primary
                r.fail_reason = fwd_fail
                r.strand = "+"
                return [r]

            # RC pass succeeded — switch working sequences so insert orientation
            # in the result is in mRNA (forward) orientation
            trimmed_seq  = rc_seq
            trimmed_qual = rc_qual
        else:
            # RC disabled in config — report forward failure and stop
            r.fail_reason = fwd_fail
            return [r]

    # Unpack the successful extraction dict into ReadResult fields
    r.strand            = extraction["strand"]
    r.linker1_start     = extraction["linker1_start"]
    r.linker1_end       = extraction["linker1_end"]
    r.linker1_edit_dist = extraction["linker1_edit_dist"]
    r.linker2_start     = extraction["linker2_start"]
    r.linker2_end       = extraction["linker2_end"]
    r.linker2_edit_dist = extraction["linker2_edit_dist"]
    r.polya_end         = extraction["polya_end"]
    r.bc1_raw           = extraction["bc1_raw"]
    r.bc2_raw           = extraction["bc2_raw"]
    r.bc3_raw           = extraction["bc3_raw"]
    r.umi_raw           = extraction["umi_raw"]
    r.bc1_corrected     = extraction["bc1_corrected"]
    r.bc2_corrected     = extraction["bc2_corrected"]
    r.bc3_corrected     = extraction["bc3_corrected"]
    r.bc1_edit_dist     = extraction["bc1_edit_dist"]
    r.bc2_edit_dist     = extraction["bc2_edit_dist"]
    r.bc3_edit_dist     = extraction["bc3_edit_dist"]
    r.bc1_window_len    = extraction["bc1_window_len"]
    r.bc2_window_len    = extraction["bc2_window_len"]
    r.bc3_window_len    = extraction["bc3_window_len"]
    r.insert_seq        = extraction["insert_seq"]
    r.insert_qual       = extraction["insert_qual"]

    # ------------------------------------------------------------------
    # Stage 6: Chimeric read detection
    # The insert is now in forward-strand (mRNA) orientation in both the
    # "+" and "-" cases, so adapter/linker sequences match in the forward
    # direction for both the adapter check and structure check.
    # Chimeric reads pass barcode extraction but are routed separately.
    # ------------------------------------------------------------------
    is_chimeric, chimeric_reason = detect_chimera(r.insert_seq, cfg)
    r.is_chimeric     = is_chimeric
    r.chimeric_reason = chimeric_reason
    r.chimeric_part   = 1   # Part 1 of 1 — splitting not yet implemented
    r.chimeric_total  = 1

    r.passed = True
    return [r]


# ---------------------------------------------------------------------------
# FASTQ reader (supports .gz and plain)
# ---------------------------------------------------------------------------

def read_fastq(path: str):
    """
    Generator yielding (read_id, seq, qual) tuples from a FASTQ file.

    Supports both plain text and gzip-compressed (.gz) FASTQ files.
    Sequences are uppercased on yield so all downstream comparisons are
    case-insensitive without needing per-function normalisation.

    The read_id is taken as the first whitespace-delimited token after '@'
    on the header line, matching the behaviour of standard tools like minimap2
    and samtools.

    Args:
        path: Path to the FASTQ file (.fastq or .fastq.gz).

    Yields:
        (read_id, seq, qual) tuples.
    """
    opener = gzip.open if path.endswith(".gz") else open
    mode   = "rt"
    with opener(path, mode) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq  = f.readline().strip()
            f.readline()   # Skip '+' separator line
            qual = f.readline().strip()
            read_id = header[1:].split()[0]   # Strip '@'; take first token only
            yield read_id, seq.upper(), qual


def chunk_reads(iterable, chunk_size: int):
    """
    Split an iterable into fixed-size chunks for multiprocessing batching.

    Yields lists of up to chunk_size elements. The last chunk may be shorter
    if the total count is not a multiple of chunk_size.

    Args:
        iterable:   Any iterable (typically a read_fastq() generator).
        chunk_size: Maximum number of elements per chunk.

    Yields:
        Lists of elements from iterable.
    """
    it = iter(iterable)
    while True:
        chunk = list(islice(it, chunk_size))
        if not chunk:
            break
        yield chunk


# ---------------------------------------------------------------------------
# Main pipeline orchestrator
# ---------------------------------------------------------------------------

def run_pipeline(
    config_path: str,
    input_fastq: str,
    output_dir: str,
    n_cores: int,
    chunk_size: int = 1000,
):
    """
    Main entry point. Coordinates all pipeline stages across all reads.

    Args:
        config_path:  Path to config.yaml
        input_fastq:  Path to input FASTQ (plain or .gz)
        output_dir:   Directory for all output files
        n_cores:      Number of parallel processes
        chunk_size:   Reads per multiprocessing chunk
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- Load config and whitelists ---
    cfg    = load_config(config_path)
    wl_bc1 = load_whitelist(cfg.whitelist_bc1)
    wl_bc2 = load_whitelist(cfg.whitelist_bc2)
    wl_bc3 = load_whitelist(cfg.whitelist_bc3)

    # --- Output file paths ---
    out_passed   = os.path.join(output_dir, "assigned_reads.fastq")
    out_failed   = os.path.join(output_dir, "unassigned_reads.fastq")
    out_chimeric = os.path.join(output_dir, "chimeric_reads.fastq")
    out_summary  = os.path.join(output_dir, "qc_summary.tsv")

    # --- QC counters ---
    stats = defaultdict(int)
    stats["total"] = 0

    log.info("Starting pipeline | input: %s | cores: %d", input_fastq, n_cores)

    with open(out_passed, "w")   as f_pass, \
         open(out_failed, "w")   as f_fail, \
         open(out_chimeric, "w") as f_chim, \
         multiprocessing.Pool(processes=n_cores) as pool:

        for chunk in tqdm(chunk_reads(read_fastq(input_fastq), chunk_size),
                          desc="Processing reads", unit="chunk"):

            # Build per-read argument tuples for pool.map().
            # All shared objects (cfg, whitelists) are passed by value to each
            # worker — this is safe because they are read-only after loading.
            args_list = [
                (read_id, seq, qual, cfg, wl_bc1, wl_bc2, wl_bc3)
                for read_id, seq, qual in chunk
            ]

            # process_read() returns List[ReadResult] — always length 1 currently
            results_lists = pool.map(process_read, args_list)

            for result_list in results_lists:
                stats["total"] += 1
                for r in result_list:
                    if not r.passed:
                        # Barcode extraction failed — route to unassigned file
                        stats["failed"] += 1
                        stats[f"fail_{r.fail_reason}"] += 1
                        if cfg.write_failed:
                            f_fail.write(format_failed_fastq_record(r))
                    elif r.is_chimeric:
                        # Barcode extraction passed but chimeric signal found —
                        # route to chimeric file rather than the clean output
                        stats["chimeric"] += 1
                        stats[f"chimeric_{r.chimeric_reason}"] += 1
                        f_chim.write(format_fastq_record(r))
                    else:
                        # Clean pass — write to assigned output
                        stats["passed"] += 1
                        if r.strand == "-":
                            stats["passed_reverse_strand"] += 1
                        f_pass.write(format_fastq_record(r))

    # --- Write TSV summary ---
    write_qc_summary(stats, out_summary)

    log.info("Pipeline complete.")
    log.info("  Passed reads   : %d / %d (%.1f%%)",
             stats["passed"], stats["total"],
             100 * stats["passed"] / max(stats["total"], 1))
    log.info("  Chimeric reads : %d (%.1f%%)",
             stats.get("chimeric", 0),
             100 * stats.get("chimeric", 0) / max(stats["total"], 1))
    log.info("  Output FASTQ   : %s", out_passed)
    log.info("  Chimeric FASTQ : %s", out_chimeric)
    log.info("  QC summary     : %s", out_summary)


# ---------------------------------------------------------------------------
# QC summary writer
# ---------------------------------------------------------------------------

def write_qc_summary(stats: dict, path: str):
    """
    Write a TSV summary of pipeline QC metrics.
    Columns: metric, count, percentage_of_total
    """
    total = stats.get("total", 0)

    rows = [
        ("total_reads",               stats["total"]),
        ("passed_reads",              stats.get("passed", 0)),
        ("passed_forward_strand",     stats.get("passed", 0) - stats.get("passed_reverse_strand", 0)),
        ("passed_reverse_strand",     stats.get("passed_reverse_strand", 0)),
        ("chimeric_reads",            stats.get("chimeric", 0)),
        ("chimeric_INTERNAL_ADAPTER", stats.get("chimeric_INTERNAL_ADAPTER", 0)),
        ("chimeric_SECOND_BC",        stats.get("chimeric_SECOND_BC_STRUCTURE", 0)),
        ("chimeric_BOTH",             stats.get("chimeric_INTERNAL_ADAPTER_AND_SECOND_BC_STRUCTURE", 0)),
        ("failed_reads",              stats.get("failed", 0)),
        ("fail_READ_TOO_SHORT",       stats.get("fail_READ_TOO_SHORT_AFTER_TRIMMING", 0)),
        ("fail_LINKER1_NOT_FOUND",    stats.get("fail_LINKER1_NOT_FOUND", 0)),
        ("fail_LINKER2_NOT_FOUND",    stats.get("fail_LINKER2_NOT_FOUND", 0)),
        ("fail_BC1_CORRECTION",       stats.get("fail_BC1_CORRECTION_FAILED", 0)),
        ("fail_BC2_CORRECTION",       stats.get("fail_BC2_CORRECTION_FAILED", 0)),
        ("fail_BC3_CORRECTION",       stats.get("fail_BC3_CORRECTION_FAILED", 0)),
        ("fail_INSERT_SHORT",         sum(v for k, v in stats.items()
                                         if k.startswith("fail_INSERT_TOO_SHORT"))),
    ]

    with open(path, "w") as f:
        f.write("metric\tcount\tpct_of_total\n")
        for metric, count in rows:
            pct = 100 * count / max(total, 1)
            f.write(f"{metric}\t{count}\t{pct:.2f}\n")

    log.info("QC summary written to %s", path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ONT scRNA-seq barcode extraction pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config.yaml"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input FASTQ file (plain or .gz)"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Output directory (will be created if needed)"
    )
    parser.add_argument(
        "--cores", type=int,
        default=multiprocessing.cpu_count(),
        help=f"Number of parallel processes (default: all available = {multiprocessing.cpu_count()})"
    )
    parser.add_argument(
        "--chunk_size", type=int, default=1000,
        help="Reads per processing chunk (default: 1000)"
    )

    args = parser.parse_args()

    run_pipeline(
        config_path  = args.config,
        input_fastq  = args.input,
        output_dir   = args.output_dir,
        n_cores      = args.cores,
        chunk_size   = args.chunk_size,
    )


if __name__ == "__main__":
    main()
