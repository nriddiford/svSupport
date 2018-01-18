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
    """Tests for `microhom.py`."""

    def test_bp1_read_count(self):
        """Are the correct number of bp1 supporting reads reported?"""

        bp1_reads, bp2_reads = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp1_reads == 4)

    def test_bp2_read_count(self):
        """Are the correct number of bp1 supporting reads reported?"""

        bp1_reads, bp2_reads = findSupport(bam_in, chrom, bp1, bp2, slop)
        self.assertTrue(bp2_reads == 5)


if __name__ == '__main__':
    unittest.main()
