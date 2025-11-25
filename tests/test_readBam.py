
import unittest

from svdb.readBAM import (read_bam_file, bam_to_vcf_format, SVEvidence, extract_sv_from_cigar, extract_sv_from_split_reads,)

class TestReadBAMModule(unittest.TestCase):
    def test_read_bam_nonexistent_returns_list(self):
        res = read_bam_file("this_file_does_not_exist.bam", "sample")
        self.assertIsInstance(res, list)

    def test_bam_to_vcf_format_basic(self):
        ev = SVEvidence("1", 100, "1", 200, "DEL", ["r1"])
        v = bam_to_vcf_format([ev], "mysample")
        self.assertEqual(len(v), 1)
        chromA, posA, chromB, posB, svtype, INFO, FORMAT, sample = v[0]
        self.assertEqual(chromA, "1")
        self.assertEqual(posA, 100)
        self.assertEqual(svtype, "DEL")
        self.assertIn("SVTYPE", INFO)
        self.assertEqual(INFO["SVTYPE"], "DEL")
        self.assertIn("SVLEN", INFO)
        self.assertEqual(FORMAT["GT"], ["0/1"])

    def test_extract_sv_from_cigar_on_missing_file(self):

        res = extract_sv_from_cigar("nope.bam", min_indel_size=1000)
        self.assertEqual(res, [])

    def test_extract_sv_from_split_reads_on_missing_file(self):
        res = extract_sv_from_split_reads("nope.bam", min_mapq=60, min_sv_size=1000)
        self.assertEqual(res, [])


if __name__ == "__main__":
    unittest.main()

