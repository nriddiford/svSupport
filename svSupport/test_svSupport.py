import sys
sys.dont_write_bytecode = True

import unittest

from svSupport import *

bam_in = '../data/test.bam'
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


class opposing_reads(unittest.TestCase):
    """Test reads returned that oppose breakpoints"""

    def test_bp1_opposing_read_count(self):
        """Are the correct number of bp1 opposing reads reported?"""
        bp1_refuting_reads, bp1_refuting_read_count, bp2_refuting_reads, bp2_refuting_read_count = findOpposing(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp1_refuting_read_count == 5)


    def test_bp1_opposing_read_names(self):
        """Are the correct bp1 supporting reads reported?"""
        bp1_true_opposing = ['HWI-D00405:129:C6KNAANXX:4:1308:17331:61532', 'HWI-D00405:129:C6KNAANXX:4:2104:6715:89377', 'HWI-D00405:129:C6KNAANXX:4:1304:5668:41607', 'HWI-D00405:129:C6KNAANXX:4:2306:20091:19131', 'HWI-D00405:129:C6KNAANXX:4:1206:21081:69455']

        bp1_refuting_reads, bp1_refuting_read_count, bp2_refuting_reads, bp2_refuting_read_count = findOpposing(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(sorted(bp1_refuting_reads) == sorted(bp1_true_opposing))

    def test_bp2_opposing_read_count(self):
        """Are the correct number of bp1 opposing reads reported?"""
        bp1_refuting_reads, bp1_refuting_read_count, bp2_refuting_reads, bp2_refuting_read_count = findOpposing(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp2_refuting_read_count == 11)


    def test_bp2_opposing_read_names(self):
        """Are the correct bp1 supporting reads reported?"""
        bp2_true_opposing = ['HWI-D00405:129:C6KNAANXX:4:1313:4501:44764', 'HWI-D00405:129:C6KNAANXX:4:1313:12748:63630', 'HWI-D00405:129:C6KNAANXX:4:2101:7585:78240', 'HWI-D00405:129:C6KNAANXX:4:2306:8064:69868', 'HWI-D00405:129:C6KNAANXX:4:2301:4836:91092', 'HWI-D00405:129:C6KNAANXX:4:2214:19389:55608', 'HWI-D00405:129:C6KNAANXX:4:2114:6250:93372', 'HWI-D00405:129:C6KNAANXX:4:2116:20962:76910', 'HWI-D00405:129:C6KNAANXX:4:1213:10124:91137', 'HWI-D00405:129:C6KNAANXX:4:1208:18379:99092', 'HWI-D00405:129:C6KNAANXX:4:2108:10027:6512']

        bp1_refuting_reads, bp1_refuting_read_count, bp2_refuting_reads, bp2_refuting_read_count = findOpposing(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(sorted(bp2_refuting_reads) == sorted(bp2_true_opposing))



class allele_frequency(unittest.TestCase):
    """Test correct allele frequency is returned"""
    def test_allele_frequency(self):
        bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = findSupport(bam_in, chrom, bp1, bp2, slop)
        bp1_opposing_reads, bp1_opposing_read_count, bp2_opposing_reads, bp2_opposing_read_count = findOpposing(bam_in, chrom, bp1, bp2, slop)
        allele_frequency = calculate_allele_freq(bp1_read_count, bp2_read_count, bp1_opposing_read_count, bp2_opposing_read_count)
        self.assertTrue(allele_frequency == 0.36)




if __name__ == '__main__':
    unittest.main()
