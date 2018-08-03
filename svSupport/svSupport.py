#!/usr/bin/env python
from __future__ import division
import re, sys
from optparse import OptionParser

import pandas as pd

from find_reads import FindReads
from calculate_allele_freq import AlleleFrequency
from merge_bams import merge_bams, sort_bam
from count_reads import count_reads, region_depth
from utils import *
from guess_tpye import *


def parse_config(options):
    print("\nExtracting arguments from config file: %s" % options.config)
    base_name = (os.path.splitext(options.config)[0])
    if not options.variants_out:
        sample = base_name.split('_')[0]
        outfile = sample + '_svSupport.txt'
        options.variants_out = outfile

    df = pd.read_csv(options.config, delimiter="\t")
    df = df.where((pd.notnull(df)), None)

    for i in df.index:
        options.in_file = df.loc[i, 'bam']

        options.region = df.loc[i, 'position']
        options.purity = float(df.loc[i, 'tumour_purity'])
        options.normal_bam = df.loc[i, 'normal_bam']
        options.find_bps = True
        options.guess = df.loc[i, 'guess']

        # if df.loc[i, 'type'] in ['DEL', 'DUP']:
        if df.loc[i, 'chromosome1'] == df.loc[i, 'chromosome2']:
            chrom, bp1, bp2, allele_frequency = worker(options)
        else:
            chrom, bp1, bp2, allele_frequency = (0,0,0,0)

        df.loc[i, 'alf2'] = allele_frequency
        df.loc[i, 'bp1_c'] = bp1
        df.loc[i, 'bp2_c'] = bp2

    df = df.drop(['bam', 'normal_bam', 'tumour_purity', 'guess', 'sample'], axis=1)
    df.to_csv(outfile, sep="\t", index=False)


def worker(options):
    print " -> Running svSupport new"
    bam_in = options.in_file
    normal = options.normal_bam
    region = options.region
    out_dir = options.out_dir
    debug = options.debug
    purity = float(options.purity)
    find_bps = options.find_bps
    test = options.test
    variants_out = options.variants_out
    guess = options.guess

    chrom, bp1, bp2 = re.split(':|-', region)
    bp1 = int(bp1)
    bp2 = int(bp2)

    if debug:
        print_options(bam_in, normal, chrom, bp1, bp2, find_bps, debug, test, out_dir)

    if options.config:
        print("python svSupport.py -i %s -n %s -l %s:%s-%s -p %s -f %s -o %s -v %s") % (
        bam_in, normal, chrom, bp1, bp2, purity, find_bps, out_dir, variants_out)

    if guess:
        bp1_reads, bp1_best_guess = guess_type(chrom, bp1, 'bp1', options)
        bp1_best_guess = max(bp1_reads, key=bp1_reads.get)
        bp2_reads, bp2_best_guess = guess_type(chrom, bp2, 'bp2', options)
        bp2_best_guess = max(bp2_reads, key=bp2_reads.get)

        bp1_best_guess, bp2_best_guess, sv_type, read_sig = classify_sv(bp1_best_guess, bp2_best_guess, )
        print("SV type : %s" % sv_type)
        print("Read signature : %s" % read_sig)
    else:
        bp1_best_guess = 'F_bp1'
        bp2_best_guess = 'bp2_R'

    if normal:
        print("* Calculating allele frequency from read depth file: %s" % bam_in)
        opposing, supporting, adj_ratio = get_depth(bam_in, normal, chrom, bp1, bp2)

        af = AlleleFrequency(opposing, supporting, purity, chrom)
        allele_frequency, adj_ratio = af.read_depth_af()
        classify_cnv(chrom, adj_ratio)

        return (chrom, bp1, bp2, allele_frequency)
    else:
        if find_bps:
            bp1, bp1_count = hone_bps(bam_in, chrom, bp1, bp1_best_guess)
            bp2, bp2_count = hone_bps(bam_in, chrom, bp2, bp2_best_guess)

            print("* Bp1 adjusted to: %s [%s split reads found]") % (bp1, bp1_count)
            print("* Bp2 adjusted to: %s [%s split reads found]") % (bp2, bp2_count)

        make_dirs(out_dir)

        bp_regions, slop = get_regions(bam_in, chrom, bp1, bp2, out_dir)

        reads = FindReads(bp_regions, chrom, bp1, bp2, slop, out_dir, debug, bp1_best_guess, bp2_best_guess)
        bp1_supporting_reads, bp1_support_count, bp1_support_bam, bp1_opposing_reads, bp1_oppose_count, bp1_oppose_bam = reads.bp1_reads()
        bp2_supporting_reads, bp2_support_count, bp2_support_bam, bp2_opposing_reads, bp2_oppose_count, bp2_oppose_bam = reads.bp2_reads(
            bp1_supporting_reads, bp1_opposing_reads)

        s = '_'.join(map(str, [chrom, bp1, bp2, 'sv_support'])) + '.bam'
        o = '_'.join(map(str, [chrom, bp1, bp2, 'sv_oppose'])) + '.bam'

        support_out = os.path.join(out_dir, s)
        oppose_out = os.path.join(out_dir, o)

        merge_bams(support_out, out_dir, [bp1_support_bam, bp2_support_bam])
        merge_bams(oppose_out, out_dir, [bp1_oppose_bam, bp2_oppose_bam])

        all_su_reads = bp1_supporting_reads + bp2_supporting_reads
        total_support = len(set(all_su_reads))

        all_op_reads = bp1_opposing_reads + bp2_opposing_reads
        total_oppose = len(set(all_op_reads))

        print("* Found %s reads in support of variant" % total_support)
        print("* Found %s reads opposing variant" % total_oppose)

        # allele_frequency = calculate_allele_freq(total_support, total_oppose, purity)
        af = AlleleFrequency(total_oppose, total_support, purity, chrom)
        allele_frequency = af.read_support_af()

        return (chrom, bp1, bp2, allele_frequency)


def get_depth(bam_in, normal, chrom, bp1, bp2):
    chromosomes = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
    t_reads_by_chrom, tumour_mapped = count_reads(bam_in, chromosomes)
    t_read_count = region_depth(bam_in, chrom, bp1, bp2)

    n_reads_by_chrom, normal_mapped = count_reads(normal, chromosomes)
    n_read_count = region_depth(normal, chrom, bp1, bp2)

    mapped_ratio = tumour_mapped / normal_mapped

    if mapped_ratio < 1:
        t_corr = t_read_count
        n_corr = round((n_read_count * mapped_ratio))
    else:
        t_corr = round((t_read_count * mapped_ratio))
        n_corr = n_read_count

    adj_ratio = round((t_corr / n_corr), 2)
    print("Normalised read count ratio: %s (%s/%s)") % (adj_ratio, t_corr, n_corr)
    return (n_corr, t_corr, adj_ratio)


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

    samfile = pysam.Samfile(regions, "rb")
    dups_rem = os.path.join(out_dir, "bp_regions" + ".bam")

    with pysam.AlignmentFile(dups_rem, "wb", template=samfile) as out:
        for read in samfile.fetch():
            if read.is_duplicate:
                continue
            out.write(read)

    os.remove(os.path.join(out_dir, 'bp_regs' + '.s.bam'))
    os.remove(os.path.join(out_dir, 'bp_regs' + '.s.bam.bai'))

    sorted_bam = sort_bam(out_dir, dups_rem)

    return sorted_bam, slop


def get_args():
    parser = OptionParser()

    parser.add_option("-i",
                      "--in_file",
                      dest="in_file",
                      action="store",
                      help="A sorted .bam file containing the reads " +
                           "supporting the structural variant calls",
                      metavar="FILE")

    parser.add_option("-n",
                      "--normal_bam",
                      dest="normal_bam",
                      action="store",
                      help="A sorted .bam file for the normal sample " +
                           "used for calculating allele frequency based " +
                           "on read depth",
                      metavar="FILE")

    parser.add_option("-p",
                      "--purity",
                      dest="purity",
                      action="store",
                      help="Tumour purity e.g. 0.75 " +
                           "[Default: 1]")

    parser.add_option("-f",
                      "--find_bps",
                      dest="find_bps",
                      action="store_true",
                      help="Look for bps if position not exact " +
                           "[Default: False]")

    parser.add_option("-l",
                      "--loci",
                      dest="region",
                      action="store",
                      help="The chromosome and breakpoints for a " +
                           "structural variant in the format: " +
                           "'chrom:bp_1-bp_2'")

    parser.add_option("-o",
                      "--out_dir",
                      dest="out_dir",
                      action="store",
                      help="Directory to write output to " +
                           "[Default: '../out']")

    parser.add_option("-d",
                      "--debug",
                      dest="debug",
                      action="store_true",
                      help="Run in debug mode " +
                           "[Default: False]")

    parser.add_option("-t",
                      "--test",
                      dest="test",
                      action="store_true",
                      help="Run on test data")

    parser.add_option("-c",
                      "--config",
                      dest="config",
                      action="store",
                      help="Config file for batch processing ")

    parser.add_option("-v",
                      "--variants",
                      dest="variants_out",
                      action="store",
                      help="File to write parsed values to ")

    parser.add_option("-g",
                      "--guess",
                      dest="guess",
                      action="store_true",
                      help="Guess type of SV for read searching")

    # out_path = os.path.abspath('../out')
    parser.set_defaults(out_dir='out',
                        purity=1)

    options, args = parser.parse_args()

    if (options.in_file is None or options.region is None) and not options.test and not options.config:
        parser.print_help()
        print

    return (options, args)


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