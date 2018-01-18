import sys
sys.dont_write_bytecode = True

import unittest

from svSupport import *

bam_in = 'data/HUM-1.tagged.filt.SC.RG.bam'
chrom = '3L'
bp1 = 9892365
bp2 = 9894889
slop = 500

class breakpoint_reads(unittest.TestCase):
    """Test reads returned in support of breakpoints"""

    def test_bp1_read_count(self):
        """Are the correct number of bp1 supporting reads reported?"""
        bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp1_read_count == 4)


    def test_bp1_read_names(self):
        """Are the correct bp1 supporting reads reported?"""
        bp1_true_support = ['HWI-D00405:129:C6KNAANXX:4:1309:3990:94686', 'HWI-D00405:129:C6KNAANXX:4:1209:9222:18319', 'HWI-D00405:129:C6KNAANXX:4:2304:19694:29523', 'HWI-D00405:129:C6KNAANXX:4:1314:2618:18304']
        bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(sorted(bp1_sv_reads) == sorted(bp1_true_support))


    def test_bp2_read_count(self):
        """Are the correct number of bp2 supporting reads reported?"""
        bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp2_read_count == 5)


    def test_bp2_read_names(self):
        """Are the correct bp2 supporting reads reported?"""
        bp2_true_support = ['HWI-D00405:129:C6KNAANXX:4:1314:2618:18304', 'HWI-D00405:129:C6KNAANXX:4:2209:8835:88441', 'HWI-D00405:129:C6KNAANXX:4:2308:11448:52099', 'HWI-D00405:129:C6KNAANXX:4:1309:3990:94686', 'HWI-D00405:129:C6KNAANXX:4:1209:9222:18319']
        bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(sorted(bp2_sv_reads) == sorted(bp2_true_support))


if __name__ == '__main__':
    unittest.main()
