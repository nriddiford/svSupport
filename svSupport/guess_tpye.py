import pysam
from collections import defaultdict
import re, os


def f_bp(read, bp):
    if not read.is_proper_pair and (not read.is_reverse and read.reference_start + read.reference_length <= bp):
        return True


def bp_r(read, bp):
    if not read.is_proper_pair and (read.is_reverse and read.reference_start >= bp):
        return True


def bp_f(read, bp):
    if not read.is_proper_pair and (not read.is_reverse and read.reference_start >= bp):
        return True


def guess_type(chrom, bp, bp_number, options):
    out_file = os.path.join(options.out_dir, bp_number + "_classifying_reads" + ".bam")
    samfile = pysam.Samfile(options.in_file, "rb")

    with pysam.AlignmentFile(out_file, "wb", template=samfile) as bpReads:
        start = bp - 200
        stop = bp + 200
        count = 0
        sv_reads = defaultdict(int)
        for read in samfile.fetch(chrom, start, stop):

            if bp == read.reference_start + 1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
                    # print("Read clipped to right: %s") % (read.cigarstring)
                    if bp_number == 'bp1':
                        sv_reads['F_bp1'] += 1
                        bpReads.write(read)
                    else:
                        sv_reads['F_bp2'] += 1
                        bpReads.write(read)
                elif re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
                    # print("Read clipped to left: %s") % (read.cigarstring)
                    if bp_number == 'bp2':
                        sv_reads['bp2_R'] += 1
                        bpReads.write(read)
                    else:
                        sv_reads['bp1_R'] += 1
                        bpReads.write(read)

            elif bp_number == 'bp1':
                if f_bp(read, bp):
                    sv_reads['F_bp1'] += 1
                    bpReads.write(read)
                elif bp_r(read, bp):
                    sv_reads['bp1_R'] += 1
                    bpReads.write(read)

            elif bp_number == 'bp2':
                if bp_r(read, bp):
                    sv_reads['bp2_R'] += 1
                    bpReads.write(read)
                elif f_bp(read, bp):
                    sv_reads['F_bp2'] += 1
                    bpReads.write(read)
                elif bp_f(read, bp):
                    sv_reads['bp2_F'] += 1
                    bpReads.write(read)

    pysam.index(out_file)
    sv_reads['NA'] = 0
    maxValKey = max(sv_reads, key=sv_reads.get)

    return (sv_reads, maxValKey)
