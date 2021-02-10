import pysam
from collections import defaultdict
from trackReads import TrackReads
from getReads import filterContamination, leftClipped, rightClipped
import os


def filter_reads(bp_regions, bp1, bp2, c1, c2, sv_type, options, supporting, opposing, bp1_sig, bp2_sig, read_tags):
    clean_reads = os.path.join(options.out_dir, 'clean_disc.bam')
    disc_support = defaultdict(int)
    split_support = defaultdict(int)
    regions = pysam.Samfile(bp_regions, "rb")

    seen_read = defaultdict(int)
    duplicates = defaultdict(int)
    supplementary_clipped = []

    with pysam.AlignmentFile(clean_reads, "wb", template=regions) as cleaned:
        for read in regions.fetch():

            if read.is_supplementary:
                read_key = '_'.join([read.query_name, 'SA'])
            else:
                read_key = read.query_name

            if read_key in seen_read:
                continue
            seen_read[read_key] += 1

            read_is_clipped = False
            clipped_supporting = False

            # TODO Probably better to do this later, and remove contaminated read & mate from supporting_reads
            contaminated_bp1, dummy = filterContamination(read, bp1, options)
            contaminated_bp2, dummy = filterContamination(read, bp2, options)
            if contaminated_bp1 or contaminated_bp2:
                supporting = supporting_remove(read, supporting, options, 'Contaminated read at breakpoint')
                continue

            for tag in read_tags[read.query_name]:
                if 'opposing' in tag:
                    supporting = supporting_remove(read, supporting, options, 'read name found in opposing reads')
                    continue
                if 'clipped' in tag:
                    mate = get_mate(read, regions)
                    read_is_clipped = True
                    read.set_tag('SV', tag, value_type='Z')
                    if clipped_support(read, bp1, 'bp1', bp1_sig, options):
                        clipped_supporting = True
                    elif mate and clipped_support(mate, bp1, 'bp1', bp1_sig, options):
                        clipped_supporting = True

                    if clipped_support(read, bp2, 'bp2', bp2_sig, options):
                        clipped_supporting = True
                    elif mate and clipped_support(mate, bp2, 'bp2', bp2_sig, options):
                        clipped_supporting = True

            if read_is_clipped and not clipped_supporting and not read.query_name in supplementary_clipped:
                supporting = supporting_remove(read, supporting, options, 'clipped in wrong direction')
                continue

            if read.is_proper_pair and not read_is_clipped and not read.query_name in supplementary_clipped:
                supporting = supporting_remove(read, supporting, options, 'not discordant or clipped')
                continue
            else:
                mate = get_mate(read, regions)

            # This should be handled earlier (using tags)
            if read.query_name in opposing:
                supporting = supporting_remove(read, supporting, options, 'read name found in opposing reads')
                continue

            dupObj = TrackReads(read, mate, c1, c2, duplicates)

            if mate is None:
                duplicates, is_dup = dupObj.check_for_clipped_dup_no_mate()
                if is_dup:
                    supporting = supporting_remove(read, supporting, options, 'Clipped duplicate')
                    continue
            else:
                mapped_chrom1 = read.reference_name
                mapped_chrom2 = mate.reference_name
                duplicates, is_dup = dupObj.check_for_standard_dup()
                if is_dup:
                    supporting = supporting_remove(read, supporting, options, 'Duplicate read')
                    continue
                duplicates, is_dup = dupObj.check_for_disc_dup()
                if is_dup:
                    supporting = supporting_remove(read, supporting, options, 'Duplicate read')
                    continue
                duplicates, is_dup = dupObj.check_for_clipped_dup()
                if is_dup:
                    supporting = supporting_remove(read, supporting, options, 'Clipped duplicate')
                    continue


            # if read.reference_start < bp1 and read.reference_end > bp1:
            #     supporting = supporting_remove(read, supporting, options)
            #     print("Overlapping read")
            #     continue
            # elif read.reference_start < bp2 and read.reference_end > bp2:
            #     supporting = supporting_remove(read, supporting, options)
            #     print("Overlapping read")
            #     continue
            #
            # if mate.reference_start < bp1 and mate.reference_end > bp1:
            #     supporting = supporting_remove(mate, supporting, options)
            #     print("Overlapping mate")
            #     continue
            # if mate.reference_start < bp2 and mate.reference_end > bp2:
            #     supporting = supporting_remove(mate, supporting, options)
            #     print("Overlapping mate")
            #     continue

            if clipped_supporting:
                split_support[read.query_name] += 1
                supporting_add(read, supporting, options, 'clipped read supporting breakpoint')
            elif not mate:
                continue

            elif sv_type == 'DEL':
                if read.reference_end <= bp1 and mate.reference_start >= bp2:
                    supporting_add(read, supporting, options, 'discordant read pair supporting DEL')
                    disc_support[read.query_name] += 1
                else:
                    supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with DEL')
                    continue

            # TODO More configurations possible?
            elif sv_type == 'BND':
                if bp1_sig == 'r_bp1' and bp2_sig == 'r_bp2':
                    if read.reference_start < bp1 and not mate.is_reverse and mate.reference_start <= bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting BND')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with BND')
                        continue

                if bp1_sig == 'bp1_r' and bp2_sig == 'bp2_r':
                    if read.is_reverse and read.reference_start >= bp1 and mate.is_reverse and mate.reference_start > bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting BND')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with BND')
                        continue
            # TODO More configurations possible...
            elif sv_type == 'TANDUP':
                if bp1_sig == 'bp1_r' and bp2_sig == 'r_bp2':
                    if read.is_reverse and read.reference_start > bp1 and not mate.is_reverse and mate.reference_end <= bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting TANDUP')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with TANDUP')
                        continue

            elif sv_type == 'TRA' and mapped_chrom1 == c1 and mapped_chrom2 == c2:
                if bp1_sig == 'bp1_r' and bp2_sig == 'bp2_r':
                    if read.is_reverse and read.reference_start >= bp1 and mate.is_reverse and mate.reference_start > bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting TRA')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with TRA')
                        continue
                elif bp1_sig == 'bp1_r' and bp2_sig == 'r_bp2':
                    if read.reference_end <= bp1 and not mate.is_reverse and mate.reference_start < bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting TRA')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with TRA')
                        continue

                elif bp1_sig == 'r_bp1' and bp2_sig == 'r_bp2':
                    if read.reference_end <= bp1 and not mate.is_reverse and mate.reference_start < bp2:
                        supporting_add(read, supporting, options, 'discordant read pair supporting TRA')
                        disc_support[read.query_name] += 1
                    else:
                        supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with TRA')
                        continue

            else:
                supporting = supporting_remove(read, supporting, options, 'discordant reads not consistent with TRA')
                continue

            if read.query_name in supporting:
                if options.debug:
                    print("[>] Writing read: %s" % (read.query_name))
                cleaned.write(read)
                if mate:
                    if options.debug:
                        print("[>] Writing mate: %s" % (read.query_name))
                    cleaned.write(mate)
            # if read.query_name in supplementary_clipped:
            #     if options.debug: print("-> Writing read: %s" % (read.query_name))
            #     cleaned.write(read)
            #     if mate:
            #         if options.debug: print("-> Writing mate: %s" % (read.query_name))
            #         cleaned.write(mate)

    clean_disc = defaultdict(int)
    for r in disc_support:
        if r not in split_support:
            clean_disc[r] += 1

    clean_disc = set(list(clean_disc.keys()))
    split_support = set(list(split_support.keys()))

    return clean_reads, supporting, clean_disc, split_support


def clipped_support(read, bp, bp_number, bp_sig, options):
    if read.is_reverse:
        direction = 'r'
    else:
        direction = 'f'

    # if read is right-clipped
    if bp == read.reference_end:
        if rightClipped(read, direction, bp_number, options) == bp_sig:
            return True
    # if read is left-clipped
    elif bp == read.reference_start + 1:
        if leftClipped(read, direction, bp_number, options) == bp_sig:
            return True
    return False


def get_mate(read, samfile):
    pointer = samfile.tell() # pointer to the current position in the BAM file
    try:
        mate = samfile.mate(read)
    except ValueError:
        return None
    finally:
        samfile.seek(pointer) # Return the BAM file to the position of read1 in the pair
    return mate


def supporting_remove(read, su, options, reason):
    if read.query_name in su:
        su = set(su)
        su.remove(read.query_name)
        if options.debug:
            print("[!] Removing read: %s : %s" % (read.query_name, reason))

    return list(su)


def supporting_add(read, su, options, reason):
    su = list(su)
    if read.query_name not in su:
        su.append(read.query_name)
        if options.debug:
            print("[o] Adding new supporting read: %s : %s" % (read.query_name, reason))

    return su
