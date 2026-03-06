# ONT Single-Cell RNA-seq Barcode Extraction Pipeline

## Overview

This pipeline processes Oxford Nanopore Technology (ONT) cDNA reads from a polyA-capture single-cell RNA-seq library. Barcodes and UMIs are appended to the 3' end of the cDNA molecule and are extracted, corrected against known whitelists, and written to a FLAMES-style FASTQ file for downstream analysis.

Because ONT sequences from either strand of the cDNA, reads arrive in both orientations. The pipeline handles this transparently at every stage — adapter trimming, linker detection, barcode extraction, and output header formatting are all orientation-aware.

---

## Read Structure

For a **forward-strand** read (barcode region at the 3' end):

```
5' ─────────────────────────────────────────────────────────────────── 3'
   [cDNA insert] · [polyA] · [BC1:8bp] · [Linker1:12bp] · [BC2:8bp]
                           · [Linker2:12bp] · [BC3:8bp] · [UMI:10bp]
```

For a **reverse-strand** read (as sequenced by ONT off the opposite strand):

```
5' ──────────────────────────────────────────────────────────────────── 3'
   [RC(UMI)] · [RC(BC3)] · [RC(Linker2)] · [RC(BC2)] · [RC(Linker1)]
             · [RC(BC1)] · [RC(polyA)] · [RC(cDNA insert)]
```

The pipeline automatically detects and handles both orientations. The output insert is always written in 5'→3' mRNA orientation and barcodes are always corrected against the forward-strand whitelist.

---

## Pipeline Stages

| Stage | Name | Description |
|-------|------|-------------|
| 0 | Adapter trimming | Orientation-aware; tries both canonical and RC adapters at each end |
| 1 | Orientation detection | Forward attempt first; RC attempted if forward fails |
| 2 | Linker detection | K-mer seed-and-extend with Levenshtein verification |
| 3 | PolyA tail detection | Hard left anchor for BC1; also cleanly delimits the insert |
| 4 | Anchored barcode + UMI extraction | Dual-anchor windows; indel-tolerant alignment replaces fixed slicing |
| 5 | Whitelist correction | Levenshtein correction of BC1, BC2, BC3 against known lists |
| 6 | Chimeric read detection | Checks for internal adapter or second BC structure in insert |
| 7 | FASTQ output | FLAMES-style header with all barcode/UMI metadata |
| 8 | QC report | TSV of per-stage pass/fail counts |

---

## Stage Details

### Stage 0 — Adapter Trimming

ONT reads carry the same adapter sequences regardless of orientation, but they appear at different ends and as reverse complements depending on which strand the pore read. To handle this without knowing the read orientation in advance, the trimmer tries **both** the canonical adapter and its reverse complement at **each** end:

| Read end | Sequences tried |
|----------|----------------|
| 5' | `adapter_5prime` and `RC(adapter_3prime)` |
| 3' | `adapter_3prime` and `RC(adapter_5prime)` |

The best hit (lowest edit distance) at each end is accepted. This ensures the downstream orientation detection always receives a cleanly trimmed sequence.

### Stage 1 — Orientation Detection

The barcode extraction attempt (Stages 2–5) is encapsulated in a single function that operates on any oriented sequence. `process_read()` calls this function twice if needed:

1. **Forward pass**: run on the trimmed read as-is.
2. **RC pass** (if forward fails): run on `reverse_complement(trimmed_seq)`, with the quality string simply reversed (`qual[::-1]`).

The successful orientation is recorded as `+` or `-` in the output header. If both passes fail, the failure reason from the forward pass is reported.

### Stage 2 — Linker Detection (K-mer + Levenshtein)

Linker detection uses a two-pass strategy to balance speed and sensitivity:

**Pass 1 — K-mer seeding (O(n))**: A k-mer index of the known linker sequence is built. The search region of the read is scanned and any position where a read k-mer matches a linker k-mer is flagged as a candidate. This avoids running a full alignment over the entire read.

**Pass 2 — Levenshtein verification**: For each candidate, a padded window (linker length ± max_errors) around the estimated position is extracted and aligned against the linker using infix (HW) alignment. The best hit within the configured tolerance is accepted.

Linker1 is searched only in the expected 3' barcode region. Linker2 is searched only *after* the confirmed Linker1 end position, enforcing the correct structural order and preventing cross-matches.

### Stage 3 — PolyA Tail Detection

The polyA tail immediately upstream of BC1 serves two purposes:

1. **Left anchor for BC1**: Together with the confirmed Linker1 start, it gives a doubly-anchored BC1 window — essential for correct extraction when BC1 contains a 1bp indel.
2. **Insert delimiter**: The insert is trimmed at the start of the polyA run rather than at an estimated position, giving a clean poly-A-free cDNA insert.

The scanner searches right-to-left (from the expected BC1 start toward the cDNA body) to find the polyA run closest to BC1 and avoid matching internal A-rich regions of the cDNA. Up to `max_mismatches` non-A bases are tolerated to account for ONT homopolymer basecalling errors.

### Stage 4 — Anchored Barcode and UMI Extraction

This is the core extraction stage. Every barcode window has confirmed anchors on **both** sides, which is what allows indel tolerance without fixed-width slicing:

```
BC1  │  LEFT:  polya_end             │  RIGHT: linker1_start
BC2  │  LEFT:  linker1_end           │  RIGHT: linker2_start
BC3  │  LEFT:  linker2_end           │  RIGHT: read_end − umi_length
UMI  │  LEFT:  read_end − umi_length │  RIGHT: read_end
```

Within each window, all whitelist entries are aligned using HW (infix) mode. An 8bp barcode with a 1bp deletion appears as 7bp in the read — the alignment finds the correct whitelist entry regardless. The UMI is always anchored to the read end, making it independent of any barcode indels.

### Stage 5 — Whitelist Correction

For each of BC1, BC2, BC3:
1. **Exact match** against whitelist → accepted at distance 0 (fast path).
2. **Levenshtein scan**: all whitelist entries within `max_levenshtein` edits are found.
3. **Unique best match** → corrected to that entry with its edit distance.
4. **Ambiguous tie** (two entries at the same distance) → read failed.
5. **No match** within tolerance → read failed.

If no whitelist is provided for a barcode, the raw extracted sequence is passed through uncorrected. This allows the pipeline to run in exploratory mode without whitelists.

### Stage 6 — Chimeric Read Detection

After successful barcode extraction, the confirmed cDNA insert is inspected for evidence of a chimeric read using two independent tests:

**Test A — Internal adapter**: the ONT adapter sequence is searched within the insert using infix alignment with Levenshtein tolerance. An adapter inside the cDNA insert strongly indicates two reads were ligated end-to-end during library preparation.

**Test B — Second barcode structure**: Linker1 and Linker2 are searched inside the insert. If both are found in the correct order with enough flanking space for BC1/BC2/BC3/UMI, a complete second barcode structure exists — a structural confirmation of a chimera.

Chimeric reads are flagged and written to `chimeric_reads.fastq` for manual inspection. They are **not** included in `assigned_reads.fastq` to avoid contaminating the cell × gene matrix.

### Stage 7 — FLAMES-style FASTQ Output

Each passing read is written with a structured header encoding all barcode and QC metadata:

```
@<read_id>_#<part>_<strand>of<total>  CB:Z:...  CR:Z:...  UR:Z:...  [XC:Z:...]
```

**Name field components:**

| Component | Meaning |
|-----------|---------|
| `_#1` | Part index (always 1; reserved for future chimeric split support) |
| `_+of1` | Forward-strand read, 1 of 1 parts |
| `_-of1` | Reverse-strand read (barcode found after RC), 1 of 1 parts |

**SAM-style tags:**

| Tag | Meaning |
|-----|---------|
| `CB:Z:` | Corrected concatenated cell barcode (BC1+BC2+BC3, forward strand) |
| `CR:Z:` | Raw extracted barcodes, underscore-separated (BC1_BC2_BC3) |
| `UR:Z:` | Raw UMI sequence (forward strand) |
| `BC1_ED:i:` | Levenshtein distance for BC1 whitelist correction |
| `BC2_ED:i:` | Levenshtein distance for BC2 whitelist correction |
| `BC3_ED:i:` | Levenshtein distance for BC3 whitelist correction |
| `L1_ED:i:` | Levenshtein distance for Linker1 match |
| `L2_ED:i:` | Levenshtein distance for Linker2 match |
| `XC:Z:` | Chimeric flag (only present on chimeric reads); value is `CHIMERIC_<reason>` |

---

## Installation

```bash
pip install pyyaml tqdm edlib
```

`edlib` is optional but strongly recommended — it provides a 10–50× speed-up over the built-in pure-Python fallback via the Myers bit-vector algorithm.

---

## Configuration

Copy `config.yaml` and fill in your library-specific values. All parameters are described in detail in the config file itself. The key required fields are:

| Config key | Description |
|------------|-------------|
| `adapter.sequence_5prime` | ONT 5' adapter sequence (forward strand) |
| `adapter.sequence_3prime` | ONT 3' adapter sequence (forward strand) |
| `linker1.sequence` | Linker1 sequence (e.g. 12 bp) |
| `linker2.sequence` | Linker2 sequence (e.g. 12 bp) |
| `barcodes.whitelist_bc1/2/3` | Paths to whitelist files (one barcode per line) |
| `insert.min_length` | Minimum cDNA insert length in bp |

See [Config Reference](#config-reference) below for a full description of every parameter.

### Whitelist file format

Plain text, one barcode per line, no header:
```
AAACCCGG
TTGGCCAA
GGTTAACC
```

---

## Running

```bash
# Basic run (uses all available CPU cores)
python extractor.py \
  --config config.yaml \
  --input reads.fastq.gz \
  --output_dir results/

# Specify core count and chunk size
python extractor.py \
  --config config.yaml \
  --input reads.fastq.gz \
  --output_dir results/ \
  --cores 8 \
  --chunk_size 2000
```

---

## Output Files

| File | Description |
|------|-------------|
| `assigned_reads.fastq` | Passing reads with full barcode/UMI header tags |
| `chimeric_reads.fastq` | Reads that passed barcode extraction but were flagged as chimeric |
| `unassigned_reads.fastq` | Reads that failed barcode extraction; header contains `FAIL_REASON:` |
| `qc_summary.tsv` | Per-stage pass/fail counts and percentages |

### QC Summary

Columns: `metric | count | pct_of_total`

| Metric | Description |
|--------|-------------|
| `total_reads` | Total reads processed |
| `passed_reads` | Reads that passed all stages and are in `assigned_reads.fastq` |
| `passed_forward_strand` | Passed reads where barcode was found on the forward strand |
| `passed_reverse_strand` | Passed reads where barcode was found after reverse complementing |
| `chimeric_reads` | Reads routed to `chimeric_reads.fastq` |
| `chimeric_INTERNAL_ADAPTER` | Chimeric reads flagged by internal adapter (Test A only) |
| `chimeric_SECOND_BC` | Chimeric reads flagged by second BC structure (Test B only) |
| `chimeric_BOTH` | Chimeric reads flagged by both tests |
| `failed_reads` | Total failed reads |
| `fail_READ_TOO_SHORT` | Read too short after adapter trimming |
| `fail_LINKER1_NOT_FOUND` | Linker1 not found in either orientation |
| `fail_LINKER2_NOT_FOUND` | Linker2 not found after Linker1 |
| `fail_BC1/2/3_CORRECTION` | Barcode not in whitelist within tolerance |
| `fail_INSERT_SHORT` | cDNA insert shorter than `min_insert_length` |

---

## Running Tests

```bash
python test_pipeline.py
```

**50 tests** covering:
- Reverse complement utility (round-trip, N preservation)
- Adapter trimming — forward strand, reverse strand (RC adapters), both ends, error tolerance
- Linker finding — exact match, 1 substitution, 1 insertion, not found, search windowing
- PolyA detection — clean run, mismatch tolerance, absent
- Barcode window alignment — exact, deletion, insertion, ambiguous rejection
- Whitelist correction — exact, 1 edit, ambiguous, no whitelist
- Reverse strand handling — RC read detected, barcodes match whitelist, insert orientation, RC disabled
- Chimeric detection — internal adapter, second BC structure, both signals, checks disabled
- FLAMES header format — forward/reverse strand format, chimeric XC tag, qual length
- End-to-end — clean read, BC indels, insert too short, whitelist mismatch

---

## Tuning

| Parameter | Conservative | Balanced | Permissive |
|-----------|-------------|----------|------------|
| `linker_max_levenshtein` | 1 | 2 | 3 |
| `bc_max_levenshtein` | 0 | 1 | 2 |
| `kmer_size` | 10 | 8 | 6 |
| `adapter_max_errors` | 1 | 2 | 3 |
| `polya_min_run` | 8 | 6 | 4 |

**Diagnostic approach using `qc_summary.tsv`:**

- **High `fail_LINKER1_NOT_FOUND`**: increase `linker1_max_levenshtein` or decrease `kmer_size`. Also check that your linker sequences are correct.
- **High `fail_BC_CORRECTION`**: first check that the whitelist is correct and complete. Only relax `bc_max_levenshtein` as a last resort — setting it too high with a dense whitelist causes ambiguous corrections.
- **High `fail_READ_TOO_SHORT`**: check that adapter sequences are correct and not overtrimming. Verify `insert.min_length` is appropriate for your library.
- **High `chimeric_reads`**: may indicate a ligation artefact in the library. Review the chimeric FASTQ to assess whether the chimeric sequences are real.
- **Unexpectedly high `passed_reverse_strand`**: normal for ONT (typically 40–60% reverse strand). Very skewed ratios may indicate an adapter or orientation issue in the library.

---

## Config Reference

All parameters are documented below. The same descriptions are also present as comments in `config.yaml`.

### `adapter` — Adapter trimming (Stage 0)

| Parameter | Type | Description |
|-----------|------|-------------|
| `sequence_5prime` | string | The ONT 5' adapter in forward-strand orientation. The pipeline also searches for `RC(sequence_3prime)` at the 5' end automatically, so always supply the canonical forward-strand sequence. |
| `sequence_3prime` | string | The ONT 3' adapter in forward-strand orientation. Similarly, `RC(sequence_5prime)` is also tried at the 3' end. |
| `max_errors` | int | Maximum Levenshtein distance for an adapter match. **Default: 2.** Accounts for typical ONT basecalling errors at read ends. |

### `linker1` / `linker2` — Linker detection (Stage 2)

| Parameter | Type | Description |
|-----------|------|-------------|
| `sequence` | string | The full linker sequence (typically 12 bp). Must be the forward-strand sequence. |
| `kmer_size` | int | Length of k-mers used in the seeding pass. Smaller values find more candidate positions (higher recall, more noise); larger values are faster. **Default: 8.** Recommended range: 6–10. |
| `max_levenshtein` | int | Maximum edit distance for a linker hit to be accepted. **Default: 2.** Values above 3 risk spurious matches in the cDNA body. |

### `barcodes` — Barcode extraction and correction (Stages 4–5)

| Parameter | Type | Description |
|-----------|------|-------------|
| `bc1_length` | int | Length of BC1 in bp. **Default: 8.** Must match the library design. |
| `bc2_length` | int | Length of BC2 in bp. **Default: 8.** |
| `bc3_length` | int | Length of BC3 in bp. **Default: 8.** |
| `umi_length` | int | Length of the UMI in bp. **Default: 10.** The UMI is always the last `umi_length` bases of the read. |
| `bc1_max_levenshtein` | int | Maximum edit distance for BC1 whitelist correction. **Default: 1.** Setting to 2 risks ambiguous corrections unless the whitelist has high minimum Hamming distance (≥ 4). |
| `bc2_max_levenshtein` | int | Same as above for BC2. **Default: 1.** |
| `bc3_max_levenshtein` | int | Same as above for BC3. **Default: 1.** |
| `whitelist_bc1` | string | Path to BC1 whitelist file (one sequence per line). Set to `null` to skip correction (raw extracted sequence is used). |
| `whitelist_bc2` | string | Path to BC2 whitelist file. |
| `whitelist_bc3` | string | Path to BC3 whitelist file. |

### `insert` — Insert length filter (Stage 5)

| Parameter | Type | Description |
|-----------|------|-------------|
| `min_length` | int | Minimum cDNA insert length in bp after barcode removal. Reads with shorter inserts are discarded. **Recommended: 100–200 bp** for typical ONT cDNA libraries. |

### `polya` — PolyA tail detection (Stage 3)

| Parameter | Type | Description |
|-----------|------|-------------|
| `min_run` | int | Minimum number of consecutive A's (within the mismatch budget) to call a polyA tail. **Default: 6.** Too low risks false positives in A-rich cDNA regions; too high may miss short tails in degraded reads. |
| `max_mismatches` | int | Maximum non-A bases tolerated within the polyA run. **Default: 1.** ONT homopolymer basecalling errors commonly introduce single non-A bases within long A-runs. Set to 0 for strict mode. |

### `orientation` — Strand orientation detection (Stage 1)

| Parameter | Type | Description |
|-----------|------|-------------|
| `try_reverse_complement` | bool | If `true` (default), and the forward-strand barcode extraction fails at any stage, the pipeline takes the reverse complement of the trimmed read and retries. The successful strand is recorded as `+` or `-` in the output header. Set to `false` only if your reads are guaranteed to be in a single orientation. **Default: true.** |

### `chimeric` — Chimeric read detection (Stage 6)

| Parameter | Type | Description |
|-----------|------|-------------|
| `adapter_check` | bool | Enable Test A: search for the ONT adapter sequence inside the cDNA insert. An internal adapter indicates two reads were ligated end-to-end. **Default: true.** |
| `adapter_seq` | string | Adapter sequence to search for internally. Leave empty (`""`) to reuse `adapter.sequence_5prime` automatically. Useful if the internal adapter differs from the sequencing adapter. |
| `adapter_max_errors` | int | Levenshtein tolerance for the internal adapter search. **Default: 3.** Higher values catch more degraded adapters but increase false-positive rate. |
| `structure_check` | bool | Enable Test B: search for a complete second Linker1+Linker2 barcode structure inside the insert. This is a structural confirmation of a chimeric read independent of adapter sequence. **Default: true.** |

### `output` — Output options

| Parameter | Type | Description |
|-----------|------|-------------|
| `write_failed` | bool | If `true`, write reads that failed barcode extraction to `unassigned_reads.fastq`. The failure reason is embedded in the header as `FAIL_REASON:<reason>`. Useful for diagnosing pipeline performance. **Default: true.** |

---

## Dependencies

| Package | Required | Purpose |
|---------|----------|---------|
| `pyyaml` | Yes | Config file parsing |
| `edlib` | Recommended | Fast Levenshtein / Myers algorithm; 10–50× faster than the built-in fallback |
| `tqdm` | Optional | Progress bar during processing |

If `edlib` is not installed the pipeline falls back to a pure-Python implementation and will run correctly but more slowly, especially with large whitelists.
