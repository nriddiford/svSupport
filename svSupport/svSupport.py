#!/usr/bin/env python
from __future__ import division
import sys, re

from collections import defaultdict

from find_reads import FindReads
from calculate_allele_freq import AlleleFrequency
from merge_bams import merge_bams, sort_bam
from depthOps import get_depth
from parseConfig import *
from getArgs import get_args
from utils import *
from guessType import get_reads, getClipped


def worker(options):
    bam_in = options.in_file
    normal = options.normal_bam
    out_dir = options.out_dir
    debug = options.debug
    purity = float(options.purity)
    find_bps = options.find_bps

    chrom, bp1, bp2 = re.split(':|-', options.region)
    bp1 = int(bp1)
    bp2 = int(bp2)

    if debug:
        print_options(bam_in, normal, chrom, bp1, bp2, find_bps, debug, options.test, out_dir)

    if options.config:
        print("python svSupport.py -i %s -n %s -l %s:%s-%s -p %s -f %s -o %s -v %s") % (
        bam_in, normal, chrom, bp1, bp2, purity, find_bps, out_dir, options.variants_out)

    if normal:
        print("* Calculating allele frequency from read depth file: %s" % bam_in)
        opposing, supporting, adj_ratio = get_depth(bam_in, normal, chrom, bp1, bp2)

        af = AlleleFrequency(opposing, supporting, purity, chrom)
        allele_frequency, adj_ratio = af.read_depth_af()
        classify_cnv(chrom, adj_ratio)

        return (chrom, bp1, bp2, allele_frequency)

    bp_regions, slop = get_regions(bam_in, chrom, bp1, bp2, out_dir)

    if options.find_bps:
        print "Guessing bp"

        bp1, bp2 = find_breakpoints(bp_regions, chrom, bp1, bp2, options)

    forward_reads = defaultdict(int)
    reverse_reads = defaultdict(int)

    bp1_split_reads, bp1_best_guess, bp1_clipped_bam, bp1_disc_bam = get_reads(bp_regions, slop, 'bp1', chrom, bp1, bp2, options, forward_reads, reverse_reads)
    print bp1_best_guess, bp1_split_reads
    bp2_split_reads, bp2_best_guess, bp2_clipped_bam, bp2_disc_bam = get_reads(bp_regions, slop, 'bp2', chrom, bp2, bp1, options, forward_reads, reverse_reads)
    print bp2_best_guess, bp2_split_reads

    clio = os.path.join(out_dir, 'clipped_reads.bam')
    disco = os.path.join(out_dir, 'discordant_reads.bam')
    merge_bams(clio, out_dir, [bp1_clipped_bam, bp2_clipped_bam])
    merge_bams(disco, out_dir, [bp1_disc_bam, bp2_disc_bam])



    # bp1_guess = {}
        # for i in range(bp1-5, bp1+5):
        #     bp1_split_reads, bp1_best_guess, bp1_clipped_bam, bp1_disc_bam = get_reads(bp_regions, slop, 'bp1', chrom, i, bp2, options, forward_reads, reverse_reads)
        #     bp1_guess[i] = bp1_split_reads
        #
        # bp1 = max(bp1_guess, key=bp1_guess.get)
        # bp2_guess = {}
        # for i in range(bp2 - 5, bp2 + 5):
        #     bp2_split_reads, bp2_best_guess, bp2_clipped_bam, bp2_disc_bam = get_reads(bp_regions, slop, 'bp2', chrom, i, bp1, options, forward_reads, reverse_reads)
        #     bp2_guess[i] = bp2_split_reads
        #
        # bp2 = max(bp2_guess, key=bp2_guess.get)
        #
        # print("Breakpoint 1 adjusted to %s (%s split reads supporting)") % (bp1, bp1_guess[bp1])
        # print("Breakpoint 2 adjusted to %s (%s split reads supporting)") % (bp2, bp2_guess[bp2])
        #
        #
        # bp2_split_reads, bp2_best_guess, bp2_clipped_bam, bp2_disc_bam = get_reads(bp_regions, slop, 'bp2', chrom, bp2, bp1, options, forward_reads, reverse_reads)
        #
        #
        # print bp2_best_guess, bp2_split_reads
        #








        # bp1_reads, bp1_best_guess = guess_type2(chrom, bp1, bp2, options)
    #     bp1_best_guess = max(bp1_reads, key=bp1_reads.get)
    #     bp2_reads, bp2_best_guess = guess_type2(chrom, bp2, 'bp2', bp1, options)
    #     bp2_best_guess = max(bp2_reads, key=bp2_reads.get)
    #
    #     bp1_best_guess, bp2_best_guess, sv_type, read_sig = classify_sv(bp1_best_guess, bp2_best_guess, )
    #     print("SV type : %s" % sv_type)
    #     print("Read signature : %s" % read_sig)
    # else:
    #     bp1_best_guess = 'F_bp1'
    #     bp2_best_guess = 'bp2_R'
    #
    # if find_bps:
    #     bp1, bp1_count = hone_bps(bam_in, chrom, bp1, bp1_best_guess)
    #     bp2, bp2_count = hone_bps(bam_in, chrom, bp2, bp2_best_guess)
    #
    #     print("* Bp1 adjusted to: %s [%s split reads found]") % (bp1, bp1_count)
    #     print("* Bp2 adjusted to: %s [%s split reads found]") % (bp2, bp2_count)

    #     make_dirs(out_dir)
    #
    #     bp_regions, slop = get_regions(bam_in, chrom, bp1, bp2, out_dir)
    #
    #     reads = FindReads(bp_regions, chrom, bp1, bp2, slop, out_dir, debug, bp1_best_guess, bp2_best_guess)
    #     bp1_supporting_reads, bp1_support_count, bp1_support_bam, bp1_opposing_reads, bp1_oppose_count, bp1_oppose_bam = reads.bp1_reads()
    #     bp2_supporting_reads, bp2_support_count, bp2_support_bam, bp2_opposing_reads, bp2_oppose_count, bp2_oppose_bam = reads.bp2_reads(
    #         bp1_supporting_reads, bp1_opposing_reads)
    #
    #     s = '_'.join(map(str, [chrom, bp1, bp2, 'sv_support'])) + '.bam'
    #     o = '_'.join(map(str, [chrom, bp1, bp2, 'sv_oppose'])) + '.bam'
    #
    #     support_out = os.path.join(out_dir, s)
    #     oppose_out = os.path.join(out_dir, o)
    #
    #     merge_bams(support_out, out_dir, [bp1_support_bam, bp2_support_bam])
    #     merge_bams(oppose_out, out_dir, [bp1_oppose_bam, bp2_oppose_bam])
    #
    #     all_su_reads = bp1_supporting_reads + bp2_supporting_reads
    #     total_support = len(set(all_su_reads))
    #
    #     all_op_reads = bp1_opposing_reads + bp2_opposing_reads
    #     total_oppose = len(set(all_op_reads))
    #
    #     print("* Found %s reads in support of variant" % total_support)
    #     print("* Found %s reads opposing variant" % total_oppose)
    #
    #     # allele_frequency = calculate_allele_freq(total_support, total_oppose, purity)
    #     af = AlleleFrequency(total_oppose, total_support, purity, chrom)
    #     allele_frequency = af.read_support_af()
    #
    #     return (chrom, bp1, bp2, allele_frequency)



def find_breakpoints(regions, chrom, bp1, bp2, options):
    samfile = pysam.Samfile(regions, "rb")
    bp1_guess = {}
    for i in range(bp1 - 5, bp1 + 5):
        split_reads = 0
        readSig = defaultdict(int)
        for read in samfile.fetch(chrom, bp1 - 10, bp1 + 10):
            readSig, split_reads, bpID = getClipped(read, i, 'f', 'bp1', readSig, split_reads, options)
        bp1_guess[i] = split_reads
    bp1 = max(bp1_guess, key=bp1_guess.get)
    print("Breakpoint 1 adjusted to %s (%s split reads supporting)") % (bp1, bp1_guess[bp1])

    bp2_guess = {}
    for i in range(bp2 - 5, bp2 + 5):
        split_reads = 0
        readSig = defaultdict(int)
        for read in samfile.fetch(chrom, bp2 - 10, bp2 + 10):
            readSig, split_reads, bpID = getClipped(read, i, 'f', 'bp2', readSig, split_reads, options)
        bp2_guess[i] = split_reads
    bp2 = max(bp2_guess, key=bp2_guess.get)
    print("Breakpoint 2 adjusted to %s (%s split reads supporting)") % (bp2, bp2_guess[bp2])


    return bp1, bp2

    # for i in range(bp1-5, bp1+5):
        #
        #
        #
        #
        #
        #
        #     bp1_split_reads, bp1_best_guess, bp1_clipped_bam, bp1_disc_bam = get_reads(bp_regions, slop, 'bp1', chrom, i, bp2, options, forward_reads, reverse_reads)
        #     bp1_guess[i] = bp1_split_reads
        #
        # bp1 = max(bp1_guess, key=bp1_guess.get)
        # bp2_guess = {}
        # for i in range(bp2 - 5, bp2 + 5):
        #     bp2_split_reads, bp2_best_guess, bp2_clipped_bam, bp2_disc_bam = get_reads(bp_regions, slop, 'bp2', chrom, i, bp1, options, forward_reads, reverse_reads)
        #     bp2_guess[i] = bp2_split_reads
        #
        # bp2 = max(bp2_guess, key=bp2_guess.get)
        #
        # print("Breakpoint 1 adjusted to %s (%s split reads supporting)") % (bp1, bp1_guess[bp1])
        # print("Breakpoint 2 adjusted to %s (%s split reads supporting)") % (bp2, bp2_guess[bp2])









def hone_bps(bam_in, chrom, bp, bp_class):
    samfile = pysam.Samfile(bam_in, "rb")
    start = bp - 5
    stop = bp + 5
    original_breakpoint = int(bp)

    bp_reads = {}
    for i in range(start, stop):
        count = 0

        for read in samfile.fetch(chrom, i - 1, i + 1):
            read_end_pos = read.reference_start + read.alen
            read_start_pos = read.reference_start

            if bp_class == 'F_bp1' or bp_class == 'F_bp2':
                if read_end_pos == i:
                    count += 1
                    bp_reads[i] = count

            elif bp_class == 'bp2_R' or bp_class == 'bp1_R':
                if read.reference_start + 1 == i:
                    count += 1
                    bp_reads[i] = count

    try:
        maxValKey = max(bp_reads, key=bp_reads.get)
        read_count = bp_reads[maxValKey]
    except ValueError:
        print("Can't find split read support - using priovided breakpoints")
        maxValKey = original_breakpoint
        read_count = 0

    return (maxValKey, read_count)


def get_regions(bam_in, chrom, bp1, bp2, out_dir):
    slop = find_is_sd(bam_in, 10000)
    samfile = pysam.Samfile(bam_in, "rb")
    bp1_bam = os.path.join(out_dir, "bp1_region" + ".bam")

    with pysam.AlignmentFile(bp1_bam, "wb", template=samfile) as bp1_region:
        for read in samfile.fetch(chrom, bp1 - slop, bp1 + slop):
            bp1_region.write(read)

    bp2_bam = os.path.join(out_dir, "bp2_region" + ".bam")
    with pysam.AlignmentFile(bp2_bam, "wb", template=samfile) as bp2_region:
        for read in samfile.fetch(chrom, bp2 - slop, bp2 + slop):
            bp2_region.write(read)

    bps_bam = os.path.join(out_dir, "bp_regs" + ".bam")
    regions = merge_bams(bps_bam, out_dir, [bp1_bam, bp2_bam])

    # samfile = pysam.Samfile(regions, "rb")
    # dups_rem = os.path.join(out_dir, "bp_regions" + ".bam")
    #
    # with pysam.AlignmentFile(dups_rem, "wb", template=samfile) as out:
    #     for read in samfile.fetch():
    #         if read.is_duplicate:
    #             continue
    #         out.write(read)
    dups_rem = rmDups(regions, "bp_regions.bam", out_dir)

    os.remove(os.path.join(out_dir, 'bp_regs' + '.s.bam'))
    os.remove(os.path.join(out_dir, 'bp_regs' + '.s.bam.bai'))

    sorted_bam = sort_bam(out_dir, dups_rem)

    return sorted_bam, slop


def rmDups(bamfile, outfile, out_dir):
    samfile = pysam.Samfile(bamfile, "rb")
    dups_rem = os.path.join(out_dir, outfile)

    forward_reads = defaultdict(int)
    reverse_reads = defaultdict(int)

    with pysam.AlignmentFile(dups_rem, "wb", template=samfile) as out:
        for read in samfile.fetch():
            if read.is_duplicate:
                continue
            elif read.is_read1:
                forward_reads[read.query_name] += 1
                if forward_reads[read.query_name] > 1:
                    continue
            elif read.is_read2:
                reverse_reads[read.query_name] += 1
                if reverse_reads[read.query_name] > 1:
                    continue

            # elif read.is_read2:
            #     reverse_reads[read.query_name] += 1
            #
            # if reverse_reads[read.query_name] > 1 or forward_reads[read.query_name] > 1:
            #     continue

            out.write(read)

    return dups_rem


def main():
    options, args = get_args()
    make_dirs(options.out_dir)

    if options.config:
        cleanup(options.out_dir)
        parse_config(options)
        sys.exit()

    elif options.test:
        print
        print("* Running in test mode...")
        print

        options.region = '3L:9892365-9894889'
        options.out_dir = 'test/test_out'
        options.in_file = 'test/data/test.bam'
        options.debug = True
        options.guess = True

    if options.in_file and options.region:
        try:
            worker(options)
            # return (chrom, bp1, bp2, allele_frequency)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return


if __name__ == "__main__":
    sys.exit(main())