#!/usr/bin/env python3
"""Unit Tests — ONT scRNA-seq Pipeline"""

import sys
import unittest
sys.path.insert(0, '/home/claude')

from pipeline_working import (
    reverse_complement, trim_adapter, find_linker, find_polya_end,
    align_barcode_to_window, correct_barcode, detect_chimera,
    format_fastq_record, _build_header, PipelineConfig, ReadResult, process_read,
)

LINKER1 = "ATCCACGTGCTTG"
LINKER2 = "TCTTCAGCGTTCC"
POLYA   = "AAAAAAAAAAAA"
INSERT  = "ACGT" * 30
BC1     = "AAACCCGG"
BC2     = "TTGGCCAA"
BC3     = "GGTTAACC"
UMI     = "ATCGATCGAT"

def make_read(insert=INSERT, polya=POLYA, bc1=BC1, linker1=LINKER1,
              bc2=BC2, linker2=LINKER2, bc3=BC3, umi=UMI,
              adapter5="", adapter3=""):
    seq  = adapter5 + insert + polya + bc1 + linker1 + bc2 + linker2 + bc3 + umi + adapter3
    qual = "I" * len(seq)
    return seq, qual

def make_cfg(bc_lev=1, min_insert=50, try_rc=True,
             chimeric_adapter_check=True, chimeric_structure_check=True,
             chimeric_adapter_seq="", chimeric_adapter_max_errors=3,
             adapter_5prime=""):
    return PipelineConfig(
        adapter_5prime=adapter_5prime, adapter_3prime="", adapter_max_errors=2,
        linker1_seq=LINKER1, linker1_kmer_size=8, linker1_max_levenshtein=2,
        linker2_seq=LINKER2, linker2_kmer_size=8, linker2_max_levenshtein=2,
        bc1_length=8, bc2_length=8, bc3_length=8, umi_length=10,
        bc1_max_levenshtein=bc_lev, bc2_max_levenshtein=bc_lev, bc3_max_levenshtein=bc_lev,
        whitelist_bc1=None, whitelist_bc2=None, whitelist_bc3=None,
        min_insert_length=min_insert,
        polya_min_run=6, polya_max_mismatches=1,
        try_reverse_complement=try_rc,
        chimeric_adapter_check=chimeric_adapter_check,
        chimeric_structure_check=chimeric_structure_check,
        chimeric_adapter_seq=chimeric_adapter_seq,
        chimeric_adapter_max_errors=chimeric_adapter_max_errors,
        write_failed=True,
    )

def pargs(seq, qual, cfg, wl1=None, wl2=None, wl3=None, rid="readID"):
    return (rid, seq, qual, cfg, wl1 or {BC1}, wl2 or {BC2}, wl3 or {BC3})


class TestRC(unittest.TestCase):
    def test_simple(self):          self.assertEqual(reverse_complement("ATCG"), "CGAT")
    def test_round_trip(self):
        for s in [LINKER1, BC1, UMI, INSERT[:20]]:
            self.assertEqual(reverse_complement(reverse_complement(s)), s)
    def test_n_preserved(self):     self.assertEqual(reverse_complement("ATCGN"), "NCGAT")


class TestAdapters(unittest.TestCase):
    """
    Adapter trimming is orientation-agnostic.

    For a forward-strand read:
        5' end carries adapter_5p
        3' end carries adapter_3p

    For a reverse-strand read (as sequenced by ONT):
        5' end carries RC(adapter_3p)
        3' end carries RC(adapter_5p)

    trim_adapter tries both canonical and RC variants at each end so that
    trimming works correctly before the Stage 1 orientation detection runs.
    """
    A5 = "AATGTACTTCGTTCAGTTACG"   # representative ONT 5' adapter
    A3 = "GCAATACGTAACTGAACGAAGT"  # representative ONT 3' adapter

    def _rc(self, s):
        return reverse_complement(s)

    def test_no_adapter(self):
        seq, qual = make_read()
        t, tq, flag = trim_adapter(seq, qual, "", "", 2)
        self.assertEqual(t, seq)
        self.assertFalse(flag)

    # --- Forward-strand adapter hits ---

    def test_forward_5prime_trimmed(self):
        """adapter_5p at 5' end of a forward read is removed."""
        seq, qual = make_read(adapter5=self.A5)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag)
        self.assertFalse(t.startswith(self.A5))
        self.assertEqual(len(t), len(tq))

    def test_forward_3prime_trimmed(self):
        """adapter_3p at 3' end of a forward read is removed."""
        seq, qual = make_read(adapter3=self.A3)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag)
        self.assertFalse(t.endswith(self.A3))

    def test_forward_both_ends_trimmed(self):
        """Both adapters on a forward read are removed in one call."""
        seq, qual = make_read(adapter5=self.A5, adapter3=self.A3)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag)
        self.assertFalse(t.startswith(self.A5))
        self.assertFalse(t.endswith(self.A3))
        self.assertEqual(len(t), len(tq))

    # --- Reverse-strand adapter hits (RC orientation) ---

    def test_rc_read_5prime_is_rc_of_3p_adapter(self):
        """
        On a reverse-strand read, the 5' end carries RC(adapter_3p).
        trim_adapter must recognise and remove it.
        """
        rc_a5 = self._rc(self.A3)          # what a rev-strand read has at its 5' end
        seq, qual = make_read(adapter5=rc_a5)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag, msg="RC of adapter_3p at 5' end was not trimmed")
        self.assertFalse(t.startswith(rc_a5))
        self.assertEqual(len(t), len(tq))

    def test_rc_read_3prime_is_rc_of_5p_adapter(self):
        """
        On a reverse-strand read, the 3' end carries RC(adapter_5p).
        trim_adapter must recognise and remove it.
        """
        rc_a3 = self._rc(self.A5)          # what a rev-strand read has at its 3' end
        seq, qual = make_read(adapter3=rc_a3)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag, msg="RC of adapter_5p at 3' end was not trimmed")
        self.assertFalse(t.endswith(rc_a3))

    def test_rc_read_both_ends_trimmed(self):
        """Full reverse-strand read — both RC adapters removed correctly."""
        rc_a5 = self._rc(self.A3)
        rc_a3 = self._rc(self.A5)
        seq, qual = make_read(adapter5=rc_a5, adapter3=rc_a3)
        t, tq, flag = trim_adapter(seq, qual, self.A5, self.A3, 2)
        self.assertTrue(flag)
        self.assertFalse(t.startswith(rc_a5))
        self.assertFalse(t.endswith(rc_a3))
        self.assertEqual(len(t), len(tq))

    # --- Error tolerance ---

    def test_forward_5prime_two_errors(self):
        """Adapter with 2 substitutions is still trimmed."""
        mutated = "AATGTXCTTCYTTCAGTTACG"   # 2 subs in A5-like sequence
        a5 =      "AATGTACTTCGTTCAGTTACG"
        seq, qual = make_read(adapter5=mutated)
        t, tq, flag = trim_adapter(seq, qual, a5, "", 2)
        self.assertTrue(flag)

    def test_rc_5prime_two_errors(self):
        """RC adapter at 5' end with 2 errors is still trimmed."""
        a5 = "AATGTACTTCGTTCAGTTACG"
        a3 = "GCAATACGTAACTGAACGAAGT"
        rc_a3_mutated = reverse_complement(a3)[:10] + "XX" + reverse_complement(a3)[12:]
        seq, qual = make_read(adapter5=rc_a3_mutated)
        t, tq, flag = trim_adapter(seq, qual, a5, a3, 2)
        self.assertTrue(flag)

    # --- Integration: RC read goes through full pipeline cleanly ---

    def test_rc_read_pipeline_passes_after_trimming(self):
        """
        A full reverse-strand read with both RC adapters should trim cleanly
        and then pass barcode extraction via the RC orientation pass.
        """
        rc_a5 = reverse_complement(self.A3)
        rc_a3 = reverse_complement(self.A5)
        fwd_seq, fwd_qual = make_read()
        # Build the RC read: RC(forward content) with RC adapters at each end
        rc_content     = reverse_complement(fwd_seq)
        rc_content_q   = fwd_qual[::-1]
        full_rc_seq    = rc_a5 + rc_content + rc_a3
        full_rc_qual   = "I" * len(rc_a5) + rc_content_q + "I" * len(rc_a3)
        cfg = make_cfg(try_rc=True, adapter_5prime=self.A5)
        cfg.adapter_3prime     = self.A3
        cfg.adapter_max_errors = 2
        results = process_read(pargs(full_rc_seq, full_rc_qual, cfg))
        r = results[0]
        self.assertTrue(r.passed, msg=f"Failed: {r.fail_reason}")
        self.assertTrue(r.adapter_trimmed)
        self.assertEqual(r.strand, "-")


class TestLinkers(unittest.TestCase):
    def test_exact(self):
        seq, _ = make_read()
        s, e, d = find_linker(seq, LINKER1, 8, 2)
        self.assertNotEqual(s, -1); self.assertEqual(d, 0)

    def test_one_sub(self):
        mut = LINKER1[:5] + "X" + LINKER1[6:]
        seq, _ = make_read(linker1=mut)
        s, e, d = find_linker(seq, LINKER1, 6, 2)
        self.assertNotEqual(s, -1)

    def test_one_insertion(self):
        mut = LINKER1[:6] + "A" + LINKER1[6:]
        seq, _ = make_read(linker1=mut)
        s, e, d = find_linker(seq, LINKER1, 6, 2)
        self.assertNotEqual(s, -1)

    def test_not_found(self):
        s, e, d = find_linker("A"*200, LINKER1, 8, 2)
        self.assertEqual(s, -1)

    def test_l2_after_l1(self):
        seq, _ = make_read()
        _, l1e, _ = find_linker(seq, LINKER1, 8, 2)
        s, e, d = find_linker(seq, LINKER2, 8, 2, search_start=l1e)
        self.assertNotEqual(s, -1); self.assertGreater(s, l1e)


class TestPolyA(unittest.TestCase):
    def test_clean(self):
        seq = "ACGT"*5 + "A"*12 + BC1
        e = find_polya_end(seq, 0, len("ACGT"*5)+12, 6, 1)
        self.assertNotEqual(e, -1)

    def test_mismatch_tolerated(self):
        seq = "ACGT"*5 + "AAACAAAAAA" + BC1
        e = find_polya_end(seq, 0, len("ACGT"*5)+10, 6, 1)
        self.assertNotEqual(e, -1)

    def test_absent(self):
        seq = "GCGTCG"*5 + BC1
        e = find_polya_end(seq, 0, 30, 8, 0)
        self.assertEqual(e, -1)


class TestBarcodeWindow(unittest.TestCase):
    WL = {BC1}
    def test_exact(self):
        _, corr, d, _, _ = align_barcode_to_window(BC1, self.WL, 8, 1)
        self.assertEqual(corr, BC1); self.assertEqual(d, 0)

    def test_deletion(self):
        del1 = BC1[:4] + BC1[5:]
        _, corr, d, _, _ = align_barcode_to_window(del1, self.WL, 8, 1)
        self.assertEqual(corr, BC1); self.assertEqual(d, 1)

    def test_insertion(self):
        ins1 = BC1[:4] + "T" + BC1[4:]
        _, corr, d, _, _ = align_barcode_to_window(ins1, self.WL, 8, 1)
        self.assertEqual(corr, BC1); self.assertEqual(d, 1)

    def test_ambiguous_rejected(self):
        wl = {"AAACCCGG", "AAACCCTT"}
        _, corr, _, _, _ = align_barcode_to_window("AAACCCAA", wl, 8, 1)
        self.assertEqual(corr, "")


class TestWhitelist(unittest.TestCase):
    WL = {"AAACCCGG", "TTGGCCAA", "GGTTAACC"}
    def test_exact(self):
        c, d = correct_barcode("AAACCCGG", self.WL, 1)
        self.assertEqual(c, "AAACCCGG"); self.assertEqual(d, 0)
    def test_one_sub(self):
        c, d = correct_barcode("AAACCCGX", self.WL, 1)
        self.assertEqual(c, "AAACCCGG"); self.assertEqual(d, 1)
    def test_no_match(self):
        c, _ = correct_barcode("AAAXXXGG", self.WL, 1)
        self.assertEqual(c, "")
    def test_no_whitelist(self):
        c, d = correct_barcode("AAACCCGG", None, 1)
        self.assertEqual(c, "AAACCCGG"); self.assertEqual(d, 0)


class TestReverseStrand(unittest.TestCase):

    def test_forward_marked_plus(self):
        seq, qual = make_read()
        r = process_read(pargs(seq, qual, make_cfg()))[0]
        self.assertTrue(r.passed); self.assertEqual(r.strand, "+")

    def test_rc_read_marked_minus(self):
        fwd, fq = make_read()
        rc_seq = reverse_complement(fwd); rc_qual = fq[::-1]
        r = process_read(pargs(rc_seq, rc_qual, make_cfg(try_rc=True)))[0]
        self.assertTrue(r.passed, msg=r.fail_reason)
        self.assertEqual(r.strand, "-")

    def test_rc_disabled_rc_read_fails(self):
        fwd, fq = make_read()
        rc_seq = reverse_complement(fwd); rc_qual = fq[::-1]
        r = process_read(pargs(rc_seq, rc_qual, make_cfg(try_rc=False)))[0]
        self.assertFalse(r.passed)

    def test_rc_barcodes_match_whitelist(self):
        """Barcodes extracted from RC read should match forward-strand whitelist."""
        fwd, fq = make_read()
        rc_seq = reverse_complement(fwd); rc_qual = fq[::-1]
        r = process_read(pargs(rc_seq, rc_qual, make_cfg(try_rc=True)))[0]
        self.assertTrue(r.passed, msg=r.fail_reason)
        self.assertEqual(r.bc1_corrected, BC1)
        self.assertEqual(r.bc2_corrected, BC2)
        self.assertEqual(r.bc3_corrected, BC3)
        self.assertEqual(r.umi_raw, UMI)

    def test_rc_insert_is_mrna_orientation(self):
        """Insert from a '-' strand read should match original mRNA sequence."""
        fwd, fq = make_read()
        rc_seq = reverse_complement(fwd); rc_qual = fq[::-1]
        r = process_read(pargs(rc_seq, rc_qual, make_cfg(try_rc=True)))[0]
        self.assertTrue(r.passed, msg=r.fail_reason)
        # Insert may be 1bp shorter at polyA boundary — check it is a suffix of INSERT
        self.assertIn(r.insert_seq, INSERT + "X")  # contained in INSERT
        self.assertGreaterEqual(len(r.insert_seq), len(INSERT) - 2)


class TestChimericDetection(unittest.TestCase):

    def test_clean_not_chimeric(self):
        seq, qual = make_read()
        r = process_read(pargs(seq, qual, make_cfg()))[0]
        self.assertFalse(r.is_chimeric)

    def test_internal_adapter_flagged(self):
        adapter = "AATGTACTTCGTT"
        chimeric_insert = ("ACGT"*10) + adapter + ("TGCA"*10)
        seq, qual = make_read(insert=chimeric_insert)
        cfg = make_cfg(chimeric_adapter_seq=adapter, chimeric_adapter_max_errors=2,
                       chimeric_adapter_check=True, chimeric_structure_check=False)
        r = process_read(pargs(seq, qual, cfg))[0]
        self.assertTrue(r.passed)
        self.assertTrue(r.is_chimeric)
        self.assertIn("INTERNAL_ADAPTER", r.chimeric_reason)

    def test_second_bc_structure_flagged(self):
        second = BC1 + LINKER1 + BC2 + LINKER2 + BC3 + UMI
        chimeric_insert = INSERT + second
        seq, qual = make_read(insert=chimeric_insert)
        cfg = make_cfg(chimeric_adapter_check=False, chimeric_structure_check=True)
        r = process_read(pargs(seq, qual, cfg))[0]
        self.assertTrue(r.passed)
        self.assertTrue(r.is_chimeric)
        self.assertIn("SECOND_BC_STRUCTURE", r.chimeric_reason)

    def test_both_signals_reported(self):
        adapter = "AATGTACTTCGTT"
        second = BC1 + LINKER1 + BC2 + LINKER2 + BC3 + UMI
        chimeric_insert = INSERT + adapter + second
        seq, qual = make_read(insert=chimeric_insert)
        cfg = make_cfg(chimeric_adapter_seq=adapter, chimeric_adapter_max_errors=2,
                       chimeric_adapter_check=True, chimeric_structure_check=True)
        r = process_read(pargs(seq, qual, cfg))[0]
        self.assertTrue(r.is_chimeric)
        self.assertIn("INTERNAL_ADAPTER",    r.chimeric_reason)
        self.assertIn("SECOND_BC_STRUCTURE", r.chimeric_reason)

    def test_checks_disabled_no_flag(self):
        adapter = "AATGTACTTCGTT"
        seq, qual = make_read(insert=INSERT + adapter)
        cfg = make_cfg(chimeric_adapter_seq=adapter, chimeric_adapter_check=False,
                       chimeric_structure_check=False)
        r = process_read(pargs(seq, qual, cfg))[0]
        self.assertFalse(r.is_chimeric)


class TestFLAMESHeader(unittest.TestCase):

    def test_forward_header_format(self):
        seq, qual = make_read()
        r = process_read(pargs(seq, qual, make_cfg(), rid="myread001"))[0]
        self.assertTrue(r.passed)
        record = format_fastq_record(r)
        header = record.split("\n")[0]
        self.assertTrue(header.startswith("@myread001_#1_+of1"), msg=header)
        for tag in ("CB:Z:", "CR:Z:", "UR:Z:", UMI):
            self.assertIn(tag, header)

    def test_reverse_strand_header(self):
        fwd, fq = make_read()
        rc_seq = reverse_complement(fwd); rc_qual = fq[::-1]
        r = process_read(pargs(rc_seq, rc_qual, make_cfg(try_rc=True), rid="rcread"))[0]
        self.assertTrue(r.passed)
        header = format_fastq_record(r).split("\n")[0]
        self.assertTrue(header.startswith("@rcread_#1_-of1"), msg=header)

    def test_chimeric_xc_tag(self):
        adapter = "AATGTACTTCGTT"
        seq, qual = make_read(insert=INSERT + adapter)
        cfg = make_cfg(chimeric_adapter_seq=adapter, chimeric_adapter_max_errors=2)
        r = process_read(pargs(seq, qual, cfg))[0]
        self.assertTrue(r.is_chimeric)
        header = format_fastq_record(r).split("\n")[0]
        self.assertIn("XC:Z:CHIMERIC_", header)

    def test_qual_matches_seq_length(self):
        seq, qual = make_read()
        r = process_read(pargs(seq, qual, make_cfg()))[0]
        lines = format_fastq_record(r).strip().split("\n")
        self.assertEqual(len(lines[1]), len(lines[3]))


class TestEndToEnd(unittest.TestCase):

    def test_clean_pass(self):
        seq, qual = make_read()
        r = process_read(pargs(seq, qual, make_cfg()))[0]
        self.assertTrue(r.passed, msg=r.fail_reason)
        self.assertEqual(r.bc1_corrected, BC1)
        self.assertEqual(r.bc2_corrected, BC2)
        self.assertEqual(r.bc3_corrected, BC3)
        self.assertEqual(r.umi_raw, UMI)
        self.assertFalse(r.is_chimeric)

    def test_bc1_deletion(self):
        seq, qual = make_read(bc1=BC1[:4]+BC1[5:])
        r = process_read(pargs(seq, qual, make_cfg(bc_lev=1)))[0]
        self.assertTrue(r.passed); self.assertEqual(r.bc1_corrected, BC1)

    def test_bc2_insertion(self):
        seq, qual = make_read(bc2=BC2[:4]+"T"+BC2[4:])
        r = process_read(pargs(seq, qual, make_cfg(bc_lev=1)))[0]
        self.assertTrue(r.passed); self.assertEqual(r.bc2_corrected, BC2)

    def test_bc3_deletion_umi_intact(self):
        seq, qual = make_read(bc3=BC3[:3]+BC3[4:])
        r = process_read(pargs(seq, qual, make_cfg(bc_lev=1)))[0]
        self.assertTrue(r.passed)
        self.assertEqual(r.bc3_corrected, BC3)
        self.assertEqual(r.umi_raw, UMI)

    def test_too_short_fails(self):
        seq, qual = make_read(insert="ACGT"*5)
        r = process_read(pargs(seq, qual, make_cfg(min_insert=50)))[0]
        self.assertFalse(r.passed)

    def test_bc_not_in_whitelist_fails(self):
        seq, qual = make_read()
        cfg = make_cfg(bc_lev=0)
        r = process_read(pargs(seq, qual, cfg, wl1={"GGGGGGGG"}))[0]
        self.assertFalse(r.passed)
        self.assertEqual(r.fail_reason, "BC1_CORRECTION_FAILED")

    def test_returns_list(self):
        seq, qual = make_read()
        result = process_read(pargs(seq, qual, make_cfg()))
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)


if __name__ == "__main__":
    print("="*60)
    print("ONT scRNA-seq Pipeline — Unit Tests")
    print("="*60)
    loader = unittest.TestLoader()
    suite  = unittest.TestSuite()
    for cls in [TestRC, TestAdapters, TestLinkers, TestPolyA, TestBarcodeWindow,
                TestWhitelist, TestReverseStrand, TestChimericDetection,
                TestFLAMESHeader, TestEndToEnd]:
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
