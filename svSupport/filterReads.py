import pysam
from collections import defaultdict
from trackReads import TrackReads
import os

def filter_reads(bp_regions, bp1, bp2, sv_type, options, supporting, bp1_sig, bp2_sig):
    clean_reads = os.path.join(options.out_dir, 'clean_disc.bam')
    disc_support = 0
    regions = pysam.Samfile(bp_regions, "rb")

    seen_read = defaultdict(int)
    duplicates = defaultdict(int)

    with pysam.AlignmentFile(clean_reads, "wb", template=regions) as cleaned:
        for read in regions.fetch():
            if read.is_proper_pair:
                continue
            else:
                mate = get_mate(read, regions)

            if read.query_name in seen_read:
                continue
            seen_read[read.query_name] += 1

            if mate is None:
                if read.query_name in supporting:
                    supporting = supporting_remove(read, supporting, options)
                continue

            chrom1 = read.reference_name
            chrom2 = mate.reference_name

            dupObj = TrackReads(read, mate, chrom1, chrom2, duplicates)
            duplicates, is_dup = dupObj.check_for_standard_dup()

            if is_dup:
                if options.debug: print "Standard dup: %s", read.query_name, read.reference_start
                continue
            duplicates, is_dup = dupObj.check_for_disc_dup()
            if is_dup:
                if options.debug: print "Discordant dup: %s", read.query_name, read.reference_start
                continue

            duplicates, is_dup = dupObj.check_for_clipped_dup()
            if is_dup:
                if options.debug: print "Clipped dup dup: %s", read.query_name, read.reference_start
                continue

            if read.reference_start < bp1 and read.reference_end > bp1:
                supporting = supporting_remove(read, supporting, options)
                continue

            if mate.reference_start < bp2 and mate.reference_end > bp2:
                supporting = supporting_remove(read, supporting, options)
                continue

            if sv_type == 'DEL':
                if read.reference_end <= bp1 and mate.reference_start >= bp2:
                    cleaned.write(read)
                    cleaned.write(mate)
                    supporting_add(read,supporting, options)
                    disc_support += 1
                else:
                    supporting = supporting_remove(read, supporting, options)
                    continue

            elif sv_type == 'BND':
                if bp1_sig == 'r_bp1' and bp2_sig == 'r_bp2':
                    if read.reference_start < bp1 and not mate.is_reverse and mate.reference_start <= bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read,supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue

                if bp1_sig == 'bp1_r' and bp2_sig == 'bp2_r':
                    if read.is_reverse and read.reference_start >= bp1 and mate.is_reverse and mate.reference_start > bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read,supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue

            elif sv_type == 'TANDUP':
                if bp1_sig == 'bp1_r' and bp2_sig == 'r_bp2':
                    if read.is_reverse and read.reference_start > bp1 and not mate.is_reverse and mate.reference_end <= bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read,supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue

            elif sv_type == 'TRA':
                if bp1_sig == 'bp1_r' and bp2_sig == 'bp2_r':
                    if read.is_reverse and read.reference_start >= bp1 and mate.is_reverse and mate.reference_start > bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read, supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue
                elif bp1_sig == 'bp1_r' and bp2_sig == 'r_bp2':
                    if read.reference_end <= bp1 and not mate.is_reverse and mate.reference_start < bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read, supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue

                elif bp1_sig == 'r_bp1' and bp2_sig == 'r_bp2':
                    if read.reference_end <= bp1 and not mate.is_reverse and mate.reference_start < bp2:
                        cleaned.write(read)
                        cleaned.write(mate)
                        supporting_add(read, supporting, options)
                        disc_support += 1
                    else:
                        supporting = supporting_remove(read, supporting, options)
                        continue


            if read.query_name in supporting:
                cleaned.write(read)
                cleaned.write(mate)
            #
            # else:
            #     bp1_clean.write(bp1_read)
            #     bp2_clean.write(bp2_read)
            # elif sv_type == 'BND':

    # sort_bam(options.out_dir, clean_reads)
    # pysam.sort(clean_reads)
    # pysam.index(clean_reads)

    return clean_reads, supporting, disc_support


def get_mate(read, samfile):
    pointer = samfile.tell() # pointer to the current position in the BAM file
    try:
        mate = samfile.mate(read)
    except ValueError:
        return
    finally:
        samfile.seek(pointer) # Return the BAM file to the position of read1 in the pair
    return mate


def supporting_remove(read, su, options):
    if read.query_name in su:
        su.remove(read.query_name)
        if options.debug: print("Removing read: %s" % (read.query_name))

    return su


def supporting_add(read, su, options):
    if read.query_name not in su:
        su.append(read.query_name)
        if options.debug: print("Adding new supporting read: %s" % (read.query_name))

    return su