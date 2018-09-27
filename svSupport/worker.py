import sys, re, os
import pysam

from collections import defaultdict

from utils import print_options, find_is_sd, getChroms
from classifyEvent import classify_sv, classify_cnv
from getReads import get_reads
from depthOps import get_depth
from findBreakpoints import find_breakpoints
from calculate_allele_freq import AlleleFrequency

from merge_bams import merge_bams, sort_bam, rm_bams


def getCooridinates(bps):
    chrom1 = None
    try:
        chrom1, bp1, bp2 = re.split(':|-', bps)
        chrom2 = chrom1
    except ValueError:
        pass

    if not chrom1:
        try:
            chrom1, bp1, chrom2, bp2 = re.split(':|-', bps)
        except ValueError:
            pass

    bp1 = int(bp1)
    bp2 = int(bp2)

    return chrom1, bp1, chrom2, bp2


def worker(options):
    bam_in = options.in_file
    normal = options.normal_bam
    out_dir = options.out_dir
    debug = options.debug
    purity = float(options.purity)
    find_bps = options.find_bps

    chrom1, bp1, chrom2, bp2 = getCooridinates(options.region)

    if debug:
        print_options(bam_in, normal, chrom1, bp1, bp2, find_bps, debug, options.test, out_dir)

    if options.config: print options

    chroms = []
    if options.chromfile:
        chroms = getChroms(options)
        print(" - Marking SV reads that don't map to one of the following chromosomes: %s") % (chroms)

    if normal:
        # if options.find_bps:
        #     options.slop = 2000
        #     bp_regions, slop = get_regions(bam_in, chrom1, bp1, chrom2, bp2, out_dir, options)
        #
        #     bp1 = find_breakpoints(bp_regions, chrom1, chrom2, bp1, 'bp1', options, cn=True)
        #     bp2 = find_breakpoints(bp_regions, chrom2, chrom2, bp2, 'bp2', options, cn=True)

        print("* Calculating allele frequency from read depth file: %s" % bam_in)
        opposing, supporting, adj_ratio = get_depth(bam_in, normal, chrom1, bp1, bp2, chroms)
        af = AlleleFrequency(opposing, supporting, purity, chrom1)
        allele_frequency, adj_ratio = af.read_depth_af()
        cnv_type = classify_cnv(chrom1, adj_ratio)

        return bp1, bp2, allele_frequency, cnv_type, '-', None

    bp_regions, slop = get_regions(bam_in, chrom1, bp1, chrom2, bp2, out_dir, options)

    if options.find_bps:
        bp1 = find_breakpoints(bp_regions, chrom1, chrom2, bp1, 'bp1', options, cn=False)
        bp2 = find_breakpoints(bp_regions, chrom2, chrom2, bp2, 'bp2', options, cn=False)

    seen_reads = []
    supporting = []
    opposing = []

    bp1_clipped_bam, bp1_disc_bam, bp1_opposing_reads, alien_integrant1, te_tagged1, bp1_sig, seen_reads, supporting, opposing = get_reads(bp_regions, 'bp1', chrom1, chrom2, bp1, bp2, options, seen_reads, chroms, supporting, opposing)
    bp2_clipped_bam, bp2_disc_bam, bp2_opposing_reads, alien_integrant2, te_tagged2, bp2_sig, seen_reads, supporting, opposing = get_reads(bp_regions, 'bp2', chrom2, chrom1, bp2, bp1, options, seen_reads, chroms, supporting, opposing)

    total_support = len(set(supporting))
    total_oppose = len(set(opposing))

    notes = []
    if bp1_sig and bp2_sig:
        if chrom1 == chrom2:
            sv_type, configuration = classify_sv(bp1_sig, bp2_sig)
        else:
            sv_type, configuration = 'TRA', '-'
        print("* Variant classified as %s") % (sv_type)
    else:
        "Read signatire not found at one of the two breakpoints - unable to classify this variant"
        sv_type, configuration = '-', '-'
        notes.append("Missing bp sig")

    print("* Found %s reads in support of variant" % total_support)
    print("* Found %s reads opposing variant" % total_oppose)

    if total_support == 0:
        print "No support found for variant"
        notes.append("No supporting reads")
        allele_frequency = '-'
    elif total_support < 4:
        print "Only found %s reads supporting variant" % total_support
        notes.append("low read support=" + str(total_support))
        af = AlleleFrequency(total_oppose, total_support, purity, chrom1)
        allele_frequency = af.read_support_af()
    else:
        af = AlleleFrequency(total_oppose, total_support, purity, chrom1)
        allele_frequency = af.read_support_af()

    svID = '_'.join(map(str, [chrom1, bp1, chrom2, bp2]))

    suoout = os.path.join(out_dir, svID + '_supporting.bam')
    opout = os.path.join(out_dir, svID + '_opposing.bam')

    merge_bams(suoout, out_dir, [bp1_clipped_bam, bp2_clipped_bam, bp1_disc_bam, bp2_disc_bam])
    merge_bams(opout, out_dir, [bp1_opposing_reads, bp2_opposing_reads])

    alien1, te1 = assessIntegration(alien_integrant1, te_tagged1, 'bp1')
    alien2, te2 = assessIntegration(alien_integrant2, te_tagged2, 'bp2')
    notes = [', '.join(notes), alien1, te1, alien2, te2]

    if not any(notes):
        notes = None

    return bp1, bp2, allele_frequency, sv_type, configuration, notes


def assessIntegration(alien, te, bp_number):
    alien[None] = 0
    ak = max(alien, key=alien.get)

    if ak:
        print " * Found %s reads supporting integration of foreign DNA at %s from source %s" % (alien[ak], bp_number, ak)
        alien_string = '_'.join(map(str, [bp_number, ak, alien[ak]]))
    else:
        alien_string = None

    te[None] = 0
    tk = max(te, key=te.get)

    if tk:
        print " * Found %s %s-tagged reads at %s" % (te[tk], tk, bp_number)
        te_string = '_'.join(map(str, [bp_number, tk, te[tk]]))
    else:
        te_string = None

    return alien_string, te_string


def get_regions(bam_in, chrom1, bp1, chrom2, bp2, out_dir, options):
    if not options.slop:
        slop = find_is_sd(bam_in, 10000)
    else:
        slop = options.slop

    samfile = pysam.Samfile(bam_in, "rb")

    bp1_bam = os.path.join(out_dir, "bp1_region" + ".bam")

    with pysam.AlignmentFile(bp1_bam, "wb", template=samfile) as bp1_region:
        for read in samfile.fetch(chrom1, bp1 - slop, bp1 + slop):
            bp1_region.write(read)

    bp2_bam = os.path.join(out_dir, "bp2_region" + ".bam")
    with pysam.AlignmentFile(bp2_bam, "wb", template=samfile) as bp2_region:
        for read in samfile.fetch(chrom2, bp2 - slop, bp2 + slop):
            bp2_region.write(read)

    bps_bam = os.path.join(out_dir, "bp_regs" + ".bam")
    regions = merge_bams(bps_bam, out_dir, [bp1_bam, bp2_bam])

    bpID = '_'.join(map(str, [chrom1, bp1-slop, chrom2, bp2+slop]))
    dups_rem = rmDups(regions, bpID + "_regions.bam", out_dir)
    sorted_bam = sort_bam(out_dir, dups_rem)
    rm_bams([dups_rem])

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

    rm_bams([bamfile])

    return dups_rem
