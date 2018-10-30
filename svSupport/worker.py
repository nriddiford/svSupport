import re
from collections import defaultdict
from utils import *
from classifyEvent import classify_sv, classify_cnv
from getReads import get_reads
from depthOps import get_depth
from findBreakpoints import find_breakpoints
from calculate_allele_freq import AlleleFrequency
from filterReads import filter_reads
import json

from merge_bams import *


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
        chrom_dict = getChroms(options)
        chroms = list(chrom_dict.keys())

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

        return bp1, bp2, allele_frequency, cnv_type, '-', None, None, None

    bp_regions, options.slop = get_regions(bam_in, chrom1, bp1, chrom2, bp2, out_dir, options, chrom_dict)

    if options.find_bps:
        bp1, bp1_split_sig = find_breakpoints(bp_regions, chrom1, chrom2, bp1, 'bp1', options, cn=False)
        bp2, bp2_split_sig = find_breakpoints(bp_regions, chrom2, chrom2, bp2, 'bp2', options, cn=False)

    seen_reads = []
    supporting = []
    opposing = []
    notes = []

    bp1_clipped_bam, bp1_disc_bam, bp1_opposing_reads, alien_integrant1, te_tagged1, bp1_disc_sig, seen_reads, supporting, opposing, contaminated_reads, bp1_read_tags = get_reads(bp_regions, 'bp1', chrom1, chrom2, bp1, bp2, options, seen_reads, chroms, supporting, opposing)

    if chrom1 not in chroms:
        supporting = []
        opposing = []
        alien_integrant1.clear()

    if contaminated_reads >= 2:
        notes.append("Contamination at bp1=" + str(contaminated_reads))

    s1 = list(supporting)
    o1 = list(opposing)

    bp2_clipped_bam, bp2_disc_bam, bp2_opposing_reads, alien_integrant2, te_tagged2, bp2_disc_sig, seen_reads, supporting, opposing, contaminated_reads, bp2_read_tags = get_reads(bp_regions, 'bp2', chrom2, chrom1, bp2, bp1, options, seen_reads, chroms, supporting, opposing)

    if chrom2 not in chroms:
        supporting = s1
        opposing = o1
        alien_integrant2.clear()

    if contaminated_reads >= 2:
        notes.append("Contamination at bp2=" + str(contaminated_reads))

    classify = True

    if not bp1_disc_sig or not bp1_split_sig and not bp2_disc_sig or not bp2_split_sig:
        print("One breakpoint has no read support. Not able to classify variant")
        sv_type, configuration = '-', '-'
        notes.append("Missing bp sig")
        classify = False
        # return bp1, bp2, 0, sv_type, configuration, notes, 0, 0

    if classify:
        if chrom1 != chrom2:
            if bp1_disc_sig and bp2_disc_sig:
                sv_type, configuration, bp1_sig, bp2_sig = classify_sv(bp1_disc_sig, bp2_disc_sig)
                print("Classifying based on split reads : TRA [%s]" % configuration)

            elif bp1_split_sig and bp2_split_sig:
                sv_type, configuration, bp1_sig, bp2_sig = classify_sv(bp1_split_sig, bp2_split_sig)
                print("Classifying based on split reads : TRA [%s]" % configuration)
            else:
                sv_type, configuration = 'TRA', '-'
            sv_type = 'TRA'

        else:
            if bp1_split_sig and bp2_split_sig:
                sv_type, configuration, bp1_sig, bp2_sig = classify_sv(bp1_split_sig, bp2_split_sig)
                print("Classifying based on split reads : %s [%s]" % (sv_type, configuration))
            elif bp1_disc_sig and bp2_disc_sig:
                sv_type, configuration, bp1_sig, bp2_sig = classify_sv(bp1_disc_sig, bp2_disc_sig)
                print("Classifying based on discordant reads : %s [%s]" % (sv_type, configuration))
            else:
                print("Read signature not found at one of the two breakpoints - unable to classify this variant")
                sv_type, configuration = '-', '-'
                notes.append("Missing bp sig")

        print("Supporting reads before filtering: %s " % len(set(supporting)))

        read_tags = merge_two_dicts(bp1_read_tags, bp2_read_tags)
        print("Breakpoint signature : %s %s" % (bp1_sig, bp2_sig))
        clean_disc_bam, supporting, disc_support, split_support = filter_reads(bp_regions, bp1, bp2, chrom1, chrom2, sv_type, options, supporting, opposing, bp1_sig, bp2_sig, read_tags)

        split_support = len(split_support)
        disc_support = len(disc_support)

        rm_bams([bp1_disc_bam, bp2_disc_bam, bp1_clipped_bam, bp2_clipped_bam])


    else:
        disc_support = 0
        split_support = 0

    print("Variant is supported by %s split reads and %s discordant read pairs" % (split_support, disc_support))
    total_support = split_support + disc_support
    total_oppose = len(set(opposing))
    print("* Found %s reads in support of variant" % total_support)
    print("* Found %s reads opposing variant" % total_oppose)

    if total_support == 0:
        print("No support found for variant")
        notes.append("No supporting reads")
    elif total_support < 3:
        notes.append("low read support=" + str(total_support))

    af = AlleleFrequency(total_oppose, total_support, purity, chrom1)
    allele_frequency = af.read_support_af()

    svID = '_'.join(map(str, [chrom1, bp1, chrom2, bp2]))
    suout = os.path.join(out_dir, svID + '_supporting_dirty.bam')
    opout = os.path.join(out_dir, svID + '_opposing.bam')

    if classify:
        susorted = merge_bams(suout, out_dir, [clean_disc_bam])
    else:
        susorted = merge_bams(suout, out_dir, [bp1_clipped_bam, bp2_clipped_bam, bp1_disc_bam, bp2_disc_bam])

    merge_bams(opout, out_dir, [bp1_opposing_reads, bp2_opposing_reads])
    snodups = os.path.join(svID + '_supporting.s.bam')
    rmDups(susorted, snodups, out_dir)

    alien1, te1 = assessIntegration(alien_integrant1, te_tagged1, 'bp1')
    alien2, te2 = assessIntegration(alien_integrant2, te_tagged2, 'bp2')
    notes = [', '.join(notes), alien1, te1, alien2, te2]

    if not any(notes):
        notes = []

    return bp1, bp2, allele_frequency, sv_type, configuration, notes, split_support, disc_support


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

    if tk and te[tk] > 1:
        print " * Found %s %s-tagged reads at %s" % (te[tk], tk, bp_number)
        te_string = '_'.join(map(str, [bp_number, tk, te[tk]]))
    else:
        te_string = None

    return alien_string, te_string


def get_regions(bam_in, chrom1, bp1, chrom2, bp2, out_dir, options, chrom_dict):
    if not options.slop:
        slop = find_is_sd(bam_in, 10000)
    else:
        slop = options.slop

    downsample = False
    samfile = pysam.Samfile(bam_in, "rb")
    overlapping_windows = False

    if chrom1 == chrom2 and bp2 - slop <= bp1 + slop:
        overlapping_windows = True
        print("Windows overlap")
        bp1_window_start, bp1_window_end = check_windows(bp1, chrom1, slop, samfile, chrom_dict)
        bp2_window_start, bp2_window_end = check_windows(bp2, chrom2, slop, samfile, chrom_dict)
        bp1_window_end = bp2_window_end
        print("Regions = %s:%s-%s" %(chrom1, bp1_window_start, bp1_window_end))
    else:
        bp1_window_start, bp1_window_end = check_windows(bp1, chrom1, slop, samfile, chrom_dict)

    bp1_bam = os.path.join(out_dir, "bp1_region" + ".bam")
    if chrom1 not in chrom_dict:
        print("%s not in %s. Will not look for breakpoints in this region" % (chrom1, chrom_dict.keys()))
        downsample = True
        r_count = 0

    with pysam.AlignmentFile(bp1_bam, "wb", template=samfile) as bp1_region:
        for read in samfile.fetch(chrom1, bp1_window_start, bp1_window_end):
            if downsample:
                r_count += 1
                if r_count > 100: break
            bp1_region.write(read)

    downsample = False

    if overlapping_windows:
        bps_bam = os.path.join(out_dir, "bp_regs" + ".bam")
        regions = merge_bams(bps_bam, out_dir, [bp1_bam])
        bpID = '_'.join(map(str, [chrom1, bp1_window_start, chrom1, bp1_window_end]))
        dups_rem = rmDups(regions, bpID + "_regions.bam", out_dir)
        sorted_bam = sort_bam(out_dir, dups_rem)
        rm_bams([dups_rem])
        return sorted_bam, slop
    else:
        if chrom2 not in chrom_dict:
            print("%s not in %s. Will not look for breakpoints in this region" % (chrom2, chrom_dict.keys()))
            downsample = True
            r_count = 0

        bp2_bam = os.path.join(out_dir, "bp2_region" + ".bam")

        bp2_window_start, bp2_window_end = check_windows(bp2, chrom2, slop, samfile, chrom_dict)

        with pysam.AlignmentFile(bp2_bam, "wb", template=samfile) as bp2_region:
            for read in samfile.fetch(chrom2, bp2_window_start, bp2_window_end):
                if downsample:
                    r_count += 1
                    if r_count > 100: break
                bp2_region.write(read)

        bps_bam = os.path.join(out_dir, "bp_regs" + ".bam")
        regions = merge_bams(bps_bam, out_dir, [bp1_bam, bp2_bam])

        bpID = '_'.join(map(str, [chrom1, bp1_window_start, chrom2, bp2_window_end]))
        dups_rem = rmDups(regions, bpID + "_regions.bam", out_dir)
        sorted_bam = sort_bam(out_dir, dups_rem)
        rm_bams([dups_rem])

    return sorted_bam, slop


def check_windows(bp, chrom, slop, samfile, chrom_dict):
    """Check that the windows are within range of chrom length.
       If not, either take beginning or end of chromosome"""
    bp_window_start = bp-slop
    bp_window_end = bp+slop

    try:
        samfile.fetch(chrom, bp_window_start, bp)
    except ValueError:
        print("%s out of range for chrom %s. Adjusting to 0" % (bp_window_start, chrom))
        bp_window_start = 0

    try:
        samfile.fetch(chrom, bp, bp_window_end)
    except ValueError:
        print("%s out of range for chrom %s. Adjusting to %s" % (bp_window_end, chrom, chrom_dict[chrom]))
        bp_window_end = int(chrom_dict[chrom])

    return bp_window_start, bp_window_end


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
    index_bam(dups_rem)

    return dups_rem