import pysam
from getReads import getClipped
from collections import defaultdict
from trackReads import TrackReads

def find_breakpoints(regions, chrom, chrom2, bp, bp_number, options):
    samfile = pysam.Samfile(regions, "rb")
    bp_guess = {}

    for i in range(bp - 5, bp + 5):
        split_reads = 0
        readSig = defaultdict(int)
        duplicates = defaultdict(int)

        for read in samfile.fetch(chrom, bp - 10, bp + 10):
            if not read.infer_read_length():
                """This is a problem - and will skip over reads with no mapped mate (which also have no cigar)"""
                continue
            dupObj = TrackReads(read, chrom, chrom2, duplicates)
            duplicates, is_dup = dupObj.check_for_standard_dup()

            if is_dup:
                continue
            duplicates, is_dup = dupObj.check_for_disc_dup()

            if is_dup:
                continue
            duplicates, is_dup = dupObj.check_for_clipped_dup()

            if is_dup:
                continue

            read, readSig, split_reads, bpID = getClipped(read, i, 'f', bp_number, readSig, split_reads, options)

        bp_guess[i] = split_reads
    bp_g = max(bp_guess, key=bp_guess.get)

    if bp_guess[bp_g] - bp_guess[bp] > 2 and bp_g != bp:
        bp = bp_g
        print("%s adjusted to %s (%s split reads supporting)") % (bp_number, bp_g, bp_guess[bp])

    return bp_guess, bp