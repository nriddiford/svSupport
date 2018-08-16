import pysam
from getReads import getClipped
from collections import defaultdict


def find_breakpoints(regions, chrom, bp, bp_number, options):
    samfile = pysam.Samfile(regions, "rb")
    bp_guess = {}

    for i in range(bp - 5, bp + 5):
        split_reads = 0
        readSig = defaultdict(int)
        for read in samfile.fetch(chrom, bp - 10, bp + 10):
            if not read.infer_read_length():
                """This is a problem - and will skip over reads with no mapped mate (which also have no cigar)"""
                continue
            if read.is_duplicate:
                continue
            read, readSig, split_reads, bpID = getClipped(read, i, 'f', bp_number, readSig, split_reads, options)

        bp_guess[i] = split_reads
    bp_g = max(bp_guess, key=bp_guess.get)

    if bp_guess[bp_g] > 0 and bp_g != bp:
        bp = bp_g
        print("%s adjusted to %s (%s split reads supporting)") % (bp_number, bp_g, bp_guess[bp])

    return bp_guess, bp