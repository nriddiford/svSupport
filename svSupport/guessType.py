import pysam
from collections import defaultdict
import re, os


def f_bp(read, bp):
    if not read.is_proper_pair and (not read.is_read1 and read.reference_start + read.reference_length <= bp):
        return True


def bp_r(read, bp):
    if not read.is_proper_pair and (read.is_reverse and read.reference_start >= bp):
        return True


def bp_f(read, bp):
    if not read.is_proper_pair and (not read.is_reverse and read.reference_start >= bp):
        return True

# def guess_type(chrom, bp, bp_number, options):
#     out_file = os.path.join(options.out_dir, bp_number + "_classifying_reads" + ".bam")
#     slop = find_is_sd(options.in_file, 10000)
#
#     samfile = pysam.Samfile(options.in_file, "rb")
#
#     with pysam.AlignmentFile(out_file, "wb", template=samfile) as bpReads:
#         start = bp - slop
#         stop = bp + slop
#
#         count = 0
#         sv_reads = defaultdict(int)
#         for read in samfile.fetch(chrom, start, stop):
#
#             if bp == read.reference_start + 1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
#                 if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
#                     if bp_number == 'bp1':
#                         sv_reads['F_bp1'] += 1
#                         bpReads.write(read)
#                     else:
#                         sv_reads['F_bp2'] += 1
#                         bpReads.write(read)
#                 elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
#                     if bp_number == 'bp2':
#                         sv_reads['bp2_R'] += 1
#                         bpReads.write(read)
#                     else:
#                         sv_reads['bp1_R'] += 1
#                         bpReads.write(read)
#
#             elif bp_number == 'bp1':
#                 if f_bp(read, bp):
#                     sv_reads['F_bp1'] += 1
#                     bpReads.write(read)
#                 elif bp_r(read, bp):
#                     sv_reads['bp1_R'] += 1
#                     bpReads.write(read)
#
#             elif bp_number == 'bp2':
#                 if bp_r(read, bp):
#                     sv_reads['bp2_R'] += 1
#                     bpReads.write(read)
#                 elif f_bp(read, bp):
#                     sv_reads['F_bp2'] += 1
#                     bpReads.write(read)
#                 elif bp_f(read, bp):
#                     sv_reads['bp2_F'] += 1
#                     bpReads.write(read)
#
#     pysam.index(out_file)
#     sv_reads['NA'] = 0
#     maxValKey = max(sv_reads, key=sv_reads.get)
#
#     return (sv_reads, maxValKey)


def disc_reads(read, bp2, forward_reads, reverse_reads):
    if read.is_read1:
        forward_reads[read.query_name] += 1
        if forward_reads[read.query_name] > 1:
            return False

    elif read.is_read2:
        reverse_reads[read.query_name] += 1
        if reverse_reads[read.query_name] > 1:
            return False

    if not read.is_proper_pair and (abs(read.next_reference_start - bp2) <= 350):
        return True


def getMate(read):
    try:
        mate = samfile.mate(read)
        return mate
    except:
        return False


def guess_type(bp_regions, slop, bp_number, chrom, bp, bp2, options, forward_reads, reverse_reads):
    clipped_out = os.path.join(options.out_dir, bp_number + "_clipped_reads" + ".bam")
    disc_out = os.path.join(options.out_dir, bp_number + "_disc_reads" + ".bam")

    samfile = pysam.Samfile(bp_regions, "rb")

    with pysam.AlignmentFile(clipped_out, "wb", template=samfile) as bpReads, pysam.AlignmentFile(disc_out, "wb", template=samfile) as discReads:
        start = bp - slop
        stop = bp + slop

        print ("Looking for reads in %s region: %s:%s-%s" % (bp_number, chrom, start, stop))
        split_reads = 0
        direction = 'f'
        breakpoint_class = defaultdict(int)

        for read in samfile.fetch(chrom, start, stop):
            read_end = read.reference_start + read.reference_length

            if read.is_reverse:
                direction = 'r'

            # bpID = fClipped(read, bp, direction, bp_number)
            # if not bpID:
            #     bpID = rClipped(read, bp, direction, bp_number)
            #
            # breakpoint_class[bpID] += 1

            # Find forward-facing clipped reads
            if bp == read_end and re.findall(r'(\d+)[S|H]', read.cigarstring):
                bpID = fClipped(read, bp, direction, bp_number)
                if not bpID:
                    bpID = rClipped(read, bp, direction, bp_number)
                if bpID:
                    breakpoint_class[bpID] += 1
                    split_reads += 1
                    bpReads.write(read)

                # if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
                #     print("----> read clipped to right: %s ---|> %s") % (read.cigarstring, read.query_name)
                #     bpID = '_'.join([direction, bp_number])
                #     breakpoint_class[bpID] += 1
                #     split_reads += 1
                #     bpReads.write(read)
                #
                # elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
                #     print("----> read clipped to left: %s |---> %s") % (read.cigarstring, read.query_name)
                #     bpID = '_'.join([bp_number, direction])
                #     breakpoint_class[bpID] += 1
                #     split_reads += 1
                #     bpReads.write(read)

            # Find reverse-facing clipped reads
            if bp == read.reference_start + 1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
                    print("<--- read clipped to right: %s <---| %s") % (read.cigarstring, read.query_name)
                    bpID = '_'.join([direction, bp_number])
                    breakpoint_class[bpID] += 1
                    split_reads += 1
                    bpReads.write(read)

                elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
                    print("<--- read clipped to left: %s <|--- %s") % (read.cigarstring, read.query_name)
                    bpID = '_'.join([bp_number, direction])
                    breakpoint_class[bpID] += 1
                    split_reads += 1
                    bpReads.write(read)


            if disc_reads(read, bp2, forward_reads, reverse_reads):
                if read.is_reverse:
                    bpID = '_'.join([bp_number, direction])
                    breakpoint_class[bpID] += 1
                else:
                    bpID = '_'.join([direction, bp_number])
                    breakpoint_class[bpID] += 1
                discReads.write(read)

            # elif read.tid == read.rnext and not read.is_duplicate:


                # pointer = samfile.tell()
                # try:
                #     mate = samfile.mate(read)
                # except ValueError:
                #     continue
                # if disc_reads(read, mate, bp2):
                #     discReads.write(read)
                # # samfile.seek(pointer)


                    # discReads.write(read)
                # discreads.append(read.qname)
            # bpReads.write(read)
            # elif bp_number == 'bp2' and disc_reads(read, bp, bp2):
            #     discreads.append(read.qname)
            #     discReads.write(read)

            # # If not clipped reads then
            # elif bp_number == 'bp1':
            #     if f_bp(read, bp):
            #         sv_reads['F_bp1'] += 1
            #         bpReads.write(read)
            #         print read
            #     elif bp_r(read, bp):
            #         sv_reads['bp1_R'] += 1
            #         # bpReads.write(read)
            #
            # elif bp_number == 'bp2':
            #     if bp_r(read, bp):
            #         sv_reads['bp2_R'] += 1
            #         bpReads.write(read)
            #     elif f_bp(read, bp):
            #         sv_reads['F_bp2'] += 1
            #         bpReads.write(read)
            #     elif bp_f(read, bp):
            #         sv_reads['bp2_F'] += 1
            #         bpReads.write(read)

    pysam.index(clipped_out)
    pysam.index(disc_out)
    breakpoint_class[None] = 0
    maxValKey = max(breakpoint_class, key=breakpoint_class.get)

    return (split_reads, maxValKey, clipped_out, disc_out)



def fClipped(read, bp, direction, bp_number):
    read_end = read.reference_start + read.reference_length
    bpID = None

    if bp == read_end and re.findall(r'(\d+)[S|H]', read.cigarstring):
        if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
            print("----> read clipped to right: %s ---|> %s") % (read.cigarstring, read.query_name)
            bpID = '_'.join([direction, bp_number])

        elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
            print("----> read clipped to left: %s |---> %s") % (read.cigarstring, read.query_name)
            bpID = '_'.join([bp_number, direction])

    return bpID

def rClipped(read, bp, direction, bp_number):
    bpID = None
    
    if bp == read.reference_start + 1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
        if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
            print("<--- read clipped to right: %s <---| %s") % (read.cigarstring, read.query_name)
            bpID = '_'.join([direction, bp_number])


        elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
            print("<--- read clipped to left: %s <|--- %s") % (read.cigarstring, read.query_name)
            bpID = '_'.join([bp_number, direction])

    return bpID