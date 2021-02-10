import pysam
from getReads import getClipped
from collections import defaultdict
from trackReads import TrackReads

def find_breakpoints(regions, chrom, chrom2, bp, bp_number, options, cn):
    samfile = pysam.Samfile(regions, "rb")
    bp_guess = {}
    sv_type_guess = {}
    read_tags = defaultdict(list)

    if cn:
        window_size = 2000
    else:
        window_size = 10

    print("Looking for reads +/- %s bps surrounding %s" % (window_size, bp_number))

    for i in range(bp - window_size, bp + window_size):
        split_reads = 0
        readSig = defaultdict(int)
        duplicates = defaultdict(int)

        for read in samfile.fetch(chrom, i - 5, i + 5):
            # TODO - can probably get rid of this?
            if not read.infer_read_length():
                """This is a problem - and will skip over reads with no mapped mate (which also have no cigar)"""
                continue
            # TODO - can probably get rid of this?
            # duplicates, is_dup = dupObj.check_for_standard_dup()
            #
            # if is_dup:
            #     continue
            # duplicates, is_dup = dupObj.check_for_disc_dup()
            #
            dupObj = TrackReads(read, read, chrom, chrom2, duplicates)
            duplicates, is_dup = dupObj.check_for_clipped_dup_no_mate()
            if is_dup:
                continue

            if read.is_reverse:
                direction = 'r'
            else:
                direction = 'f'

            read, readSig, split_reads, bpID, read_tags = getClipped(read, i, direction, bp_number, readSig, split_reads, options, read_tags)

        bp_guess[i] = split_reads
        sv_type_guess[i] = readSig

    # t_g = max(sv_type_guess, key=sv_type_guess.get)
    bp_g = max(bp_guess, key=bp_guess.get)
    svtype = sv_type_guess[bp_g]

    # TODO - Need to revert back to CN split read searching
    if cn and bp_guess[bp_g] > 3:
        bp = bp_g
        print("%s adjusted to %s (%s split reads supporting)" % (bp_number, bp_g, bp_guess[bp]))

    elif not cn and bp_guess[bp_g] > bp_guess[bp] and bp_g != bp:
        bp = bp_g
        print("%s adjusted to %s (%s split reads supporting)" % (bp_number, bp_g, bp_guess[bp]))

    return bp, svtype
