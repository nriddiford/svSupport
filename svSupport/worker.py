import sys, re, os
import pysam

from collections import defaultdict

from utils import print_options, find_is_sd, getChroms
from classifyEvent import classify_sv, classify_cnv
from getReads import get_reads
from depthOps import get_depth
from findBreakpoints import find_breakpoints
from calculate_allele_freq import AlleleFrequency

from merge_bams import merge_bams, sort_bam


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

    chroms = []
    if options.chromfile:
        chroms = getChroms(options)
        print(" - Marking SV reads that don't map to one of the following chromosomes: %s") % (chroms)


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

        bp1_guess, bp1 = find_breakpoints(bp_regions, chrom, bp1, 'bp1', options)
        bp2_guess, bp2 = find_breakpoints(bp_regions, chrom, bp2, 'bp2', options)

    seen_reads = defaultdict(int)
    bp_sig = defaultdict(int)

    supporting = []
    opposing = []

    bp1_split_reads, bp1_best_guess, bp1_clipped_bam, bp1_disc_bam, bp1_opposing_reads, non_native_split_bp1, non_native_paired_bp1, te1, te_tagged1, bp_sig, seen_reads, supporting, opposing = get_reads(bp_regions, 'bp1', chrom, bp1, bp2, options, seen_reads, chroms, bp_sig, supporting, opposing)
    bp2_split_reads, bp2_best_guess, bp2_clipped_bam, bp2_disc_bam, bp2_opposing_reads, non_native_split_bp2, non_native_paired_bp2, te2, te_tagged2, bp_sig, seen_reads, supporting, opposing = get_reads(bp_regions, 'bp2', chrom, bp2, bp1, options, seen_reads, chroms, bp_sig, supporting, opposing)

    total_support = len(set(supporting))
    total_oppose = len(set(opposing))

    print("* Found %s reads in support of variant" % total_support)
    print("* Found %s reads opposing variant" % total_oppose)

    af = AlleleFrequency(total_oppose, total_support, purity, chrom)
    allele_frequency = af.read_support_af()

    if non_native_split_bp1 or non_native_paired_bp1:
        print " * Found reads supporting integration of foreign DNA at bp1 (%s split, %s paired)" % (non_native_split_bp1, non_native_paired_bp1)

    if te_tagged1[te1]:
        print " * %s %s-tagged reads found at bp1" % (te_tagged1[te1], te1)

    if non_native_split_bp2 or non_native_paired_bp2:
        print " * Found reads supporting integration of foreign DNA at bp2 (%s split, %s paired)" % (non_native_split_bp2, non_native_paired_bp2)

    if te_tagged2[te2]:
        print " * %s %s-tagged reads found at bp2" % (te_tagged2[te2], te2)


    # clio = os.path.join(out_dir, 'clipped_reads.bam')
    # disco = os.path.join(out_dir, 'discordant_reads.bam')
    # merge_bams(disco, out_dir, [bp1_disc_bam, bp2_disc_bam])

    svID = '_'.join(map(str, [chrom, bp1, bp2]))

    suoout = os.path.join(out_dir, svID + '_supporting.bam')
    opout = os.path.join(out_dir, svID + '_opposing_reads.bam')

    merge_bams(suoout, out_dir, [bp1_clipped_bam, bp2_clipped_bam, bp1_disc_bam, bp2_disc_bam])
    merge_bams(opout, out_dir, [bp1_opposing_reads, bp2_opposing_reads])

    print "Breakpoint signatures %s" % (bp_sig)

    sv_type = classify_sv(bp_sig)



    return (chrom, bp1, bp2, allele_frequency)


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

    bpID = '_'.join(map(str, [chrom, bp1-slop, bp2+slop]))

    dups_rem = rmDups(regions, bpID + "_regions.bam", out_dir)

    sorted_bam = sort_bam(out_dir, dups_rem)

    return sorted_bam, slop


def rmDups(bamfile, outfile, out_dir):
    """Skip if read is marked as duplicate, or we see a read with the same name and starting position
       (as will happen if we are merging close by regions) more than once"""

    samfile = pysam.Samfile(bamfile, "rb")
    dups_rem = os.path.join(out_dir, outfile)
    seen_reads = defaultdict(int)

    with pysam.AlignmentFile(dups_rem, "wb", template=samfile) as out:
        for read in samfile.fetch():
            read_key = '_'.join([read.query_name, str(read.reference_start)])
            seen_reads[read_key] += 1

            if seen_reads[read_key] > 1 or read.is_duplicate:
                continue

            out.write(read)

    try:
        os.remove(bamfile)
        os.remove(bamfile + ".bai")
    except OSError:
        print("Couldn't remove %s" % bamfile)
        pass

    return dups_rem
