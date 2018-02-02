#!/usr/bin/env python
from __future__ import division
import sys, os, re
import pysam
from optparse import OptionParser

out_dir = '../out'
debug = 0


def bp1_supporting_reads(bamFile, chrom, bp1, bp2, slop):
    samfile = pysam.Samfile(bamFile, "rb")
    start=bp1-slop
    bp1_reads = []
    out_file = os.path.join(out_dir, "bp1_sv_reads" + ".bam")
    bp1_sv_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
    count = 0
    for read in samfile.fetch(chrom, start, bp1):
        read_end_pos = read.reference_start + read.alen
        mate_end_pos = read.next_reference_start + read.alen

        if read.is_duplicate:
            print(read.qname)
            continue

        if not read.is_proper_pair and not read.is_reverse and read.next_reference_start > bp2:
            if debug:
                print("* bp1 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp1_sv_reads.write(read)
            bp1_reads.append(read.qname)
            count += 1

        elif read_end_pos == bp1:
            if debug:
                print("* bp1 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.reference_start, read_end_pos)
            bp1_sv_reads.write(read)
            bp1_reads.append(read.qname)
            count += 1

    bp1_sv_reads.close()
    pysam.index(out_file)

    return(bp1_reads, count)


def bp2_supporting_reads(bamFile, chrom, bp1, bp2, slop):
    samfile = pysam.Samfile(bamFile, "rb")
    end=bp2+slop
    bp2_reads = []
    out_file = os.path.join(out_dir, "bp2_sv_reads" + ".bam")
    bp2_sv_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)

    count = 0

    for read in samfile.fetch(chrom, bp2, end):
        read_end_pos = read.reference_start + read.alen
        mate_end_pos = read.next_reference_start + read.alen

        if read.is_duplicate:
            continue

        if not read.is_proper_pair and read.is_reverse and mate_end_pos < bp1:
            if debug:
                print("* bp2 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp2_sv_reads.write(read)
            bp2_reads.append(read.qname)
            count += 1

        elif read.reference_start +1 == bp2:
            if debug:
                print("* bp2 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.reference_start, read_end_pos)
            bp2_reads.append(read.qname)
            bp2_sv_reads.write(read)
            count += 1

    bp2_sv_reads.close()
    pysam.index(out_file)

    return(bp2_reads, count)


def bp_1_opposing_reads(bamFile, chrom, bp1, bp2, slop):
    samfile = pysam.Samfile(bamFile, "rb")
    start=bp1-slop
    bp1_reads = []
    read_names = set()
    out_file = os.path.join(out_dir, "bp1_opposing_reads" + ".bam")
    bp1_opposing_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
    count = 0
    for read in samfile.fetch(chrom, start, bp1):
        read_end_pos = read.reference_start + read.alen
        mate_end_pos = read.next_reference_start + read.alen

        if read.is_duplicate:
            continue

        if read.is_proper_pair and not read.is_reverse and read.next_reference_start > bp1 and not read.is_supplementary and not read_end_pos == bp1:
            if debug:
                print("* bp1 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp1_opposing_reads.write(read)
            bp1_reads.append(read.qname)
            read_names.update(read.qname)
            count += 1

        elif read.reference_start < bp1 and read_end_pos > bp1:
            if debug:
                print("* bp1 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp1_opposing_reads.write(read)
            bp1_reads.append(read.qname)
            read_names.update(read.qname)
            count += 1

    # print(read_names)
    bp1_opposing_reads.close()
    pysam.index(out_file)

    return(bp1_reads, count)


def bp_2_opposing_reads(bamFile, chrom, bp1, bp2, slop):
    samfile = pysam.Samfile(bamFile, "rb")
    end=bp2+slop
    bp2_reads = []
    out_file = os.path.join(out_dir, "bp2_opposing_reads" + ".bam")
    bp2_opposing_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
    count = 0
    for read in samfile.fetch(chrom, bp2, end):
        read_end_pos = read.reference_start + read.alen
        mate_end_pos = read.next_reference_start + read.alen

        if read.is_duplicate:
            continue

        if read.is_proper_pair and read.is_reverse and read.next_reference_start < bp2 and not read.is_supplementary and read.reference_start +1 != bp2 and read.next_reference_start +1 != bp2:
            if debug:
                print("* bp2 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp2_opposing_reads.write(read)
            bp2_reads.append(read.qname)
            count += 1

        elif read.reference_start < bp1 and read_end_pos > bp1:
            if debug:
                print("* bp2 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
            bp2_opposing_reads.write(read)
            bp2_reads.append(read.qname)
            count += 1

    bp2_opposing_reads.close()
    pysam.index(out_file)

    return(bp2_reads, count)

def merge_bams(out_file, bams):
    in_files = ', '.join(bams)
    print("Merging bam files %s into '%s'") % (in_files, out_file)
    merge_parameters = ['-f', out_file] + bams
    pysam.merge(*merge_parameters)
    # Remove individual bp files
    for bp_file in bams:
        os.remove(bp_file)
        os.remove(bp_file + ".bai")

    pysam.index(out_file)


def get_reads(bam_in, chrom, bp1, bp2, slop, type):
    if type == 'support':
        bp1_reads, bp1_read_count = bp1_supporting_reads(bam_in, chrom, bp1, bp2, slop)
        bp2_reads, bp2_read_count = bp2_supporting_reads(bam_in, chrom, bp1, bp2, slop)
        bam1 = os.path.join(out_dir, "bp1_sv_reads" + ".bam")
        bam2 = os.path.join(out_dir, "bp2_sv_reads" + ".bam")
        out = os.path.join(out_dir, "sv_support" + ".bam")
        total_support = bp1_read_count + bp2_read_count

        print("Found %s reads in support of variant" % total_support)
    else:
        bp1_reads, bp1_read_count = bp_1_opposing_reads(bam_in, chrom, bp1, bp2, slop)
        bp2_reads, bp2_read_count = bp_2_opposing_reads(bam_in, chrom, bp1, bp2, slop)
        bam1 = os.path.join(out_dir, "bp1_opposing_reads" + ".bam")
        bam2 = os.path.join(out_dir, "bp2_opposing_reads" + ".bam")
        out = os.path.join(out_dir, "sv_oppose" + ".bam")
        total_oppose = bp1_read_count + bp2_read_count

        print("Found %s reads opposing variant" % total_oppose)

    to_merge = [bam1, bam2]
    merge_bams(out, to_merge)
    return(bp1_reads, bp1_read_count, bp2_reads, bp2_read_count)


def calculate_allele_freq(bp1_read_count, bp2_read_count, bp1_opposing_read_count, bp2_opposing_read_count, tumour_purity):
    total_support =  bp1_read_count + bp2_read_count
    total_oppose = bp1_opposing_read_count + bp2_opposing_read_count

    print("Tumour puity set to %s" % tumour_purity)

    allele_frequency = float( total_support/(total_support+total_oppose) )

    if tumour_purity == 1:
        adj_allele_frequency = allele_frequency
    else:
        adj_allele_frequency = float( total_support/( total_support + (total_oppose**tumour_purity) ) )

    allele_frequency = "{:.2f}".format(allele_frequency)
    adj_allele_frequency = "{:.2f}".format(adj_allele_frequency)
    print("Adjusted allele frequency from %s to %s") % (allele_frequency, adj_allele_frequency)

    return(allele_frequency)


def make_dirs(bam_file, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def print_options(bam_in, chrom, bp1, bp2, slop, find_bps, debug, out_dir):
    options = ['Bam file', 'Chrom', 'bp1', 'bp2', 'slop', 'search_bps', 'debug', 'Out dir']
    args = [bam_in, chrom, bp1, bp2, slop, find_bps, debug, out_dir]
    print("----")
    for index, (value1, value2) in enumerate(zip(options, args)):
         print("%s: %s") % (value1, value2)
    print("----")

def search_bps(bamFile, chrom, bp, bp_number):
    samfile = pysam.Samfile(bamFile, "rb")
    start = bp - 10
    stop = bp + 10

    print("Searching for %s in region: %s:%s-%s") % (bp_number, chrom, start, stop)

    bp_reads = {}
    for i in range(start, stop):
        count = 0

        for read in samfile.fetch(chrom, i-1, i+1):
            read_end_pos   = read.reference_start + read.alen
            read_start_pos = read.reference_start

            if bp_number == 'bp1':
                if read_end_pos == i:
                    count+=1
                    bp_reads[i] = count

            else:
                if read.reference_start +1 == i:
                    count+=1
                    bp_reads[i] = count


    maxValKey = max(bp_reads, key=bp_reads.get)
    read_count = bp_reads[maxValKey]
    return(maxValKey, read_count)


def get_args():
    parser = OptionParser()

    parser.add_option("-i", \
                      "--in_file", \
                      dest="in_file",
                      action="store",
                      help="A sorted .bam file containing the reads " + \
                         "supporting the structural variant calls", \
                         metavar="FILE")

    parser.add_option("-s", \
                    "--slop", \
                    dest="slop",
                    action="store",
                    type="int",
                    help="Distance from breakpoint to look for reads " + \
                         "[Default: 500]")

    parser.add_option("-p", \
                     "--purity", \
                     dest="purity",
                     action="store",
                     help="Tumour purity e.g. 0.75 " + \
                          "[Default: 1]")

    parser.add_option("-f", \
                     "--find_bps", \
                     dest="find_bps",
                     action="store_true",
                     help="Look for bps if position not exact " + \
                          "[Default: F ]")

    parser.add_option("-l", \
                    "--loci", \
                    dest="region",
                    action="store",
                    help="The chromosome and breakpoints for a " + \
                         "structural variant in the format: " + \
                         "'chrom:bp_1-bp_2'")

    parser.add_option("-o", \
                    "--out_dir", \
                    dest="out_dir",
                    action="store",
                    help="Directory to write output to " + \
                         "[Default: '../out']")


    parser.add_option("-d", \
                    "--debug", \
                    dest="debug",
                    action="store_true",
                    help="Run in debug mode")

    parser.set_defaults(slop=500, out_dir='../out', debug=0, purity=1, find_bps=0)
    options, args = parser.parse_args()

    if options.in_file is None or options.region is None:
      parser.print_help()
      print

    return(options, args)


def main():
    options, args = get_args()

    global out_dir
    global debug

    if options.in_file is not None:
        try:
            bam_in = options.in_file
            region = options.region
            slop   = options.slop
            out_dir = options.out_dir
            debug = options.debug
            purity = float(options.purity)
            find_bps = options.find_bps

            chrom, bp1, bp2 = re.split(':|-', region)
            bp1 = int(bp1)
            bp2 = int(bp2)
            slop = int(slop)

            if debug:
                print_options(bam_in, chrom, bp1, bp2, slop, find_bps, debug, out_dir)

            if find_bps:
                bp1, bp1_count = search_bps(bam_in, chrom, bp1, 'bp1')
                bp2, bp2_count = search_bps(bam_in, chrom, bp2, 'bp2')
                print("Bp1 adjusted to: %s [%s split reads found]") % (bp1, bp1_count)
                print("Bp1 adjusted to: %s [%s split reads found]") % (bp2, bp2_count)

            make_dirs(bam_in, out_dir)
            bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = get_reads(bam_in, chrom, bp1, bp2, slop, 'support')
            bp1_opposing_reads, bp1_opposing_read_count, bp2_opposing_reads, bp2_opposing_read_count = get_reads(bam_in, chrom, bp1, bp2, slop, 'oppose')

            allele_frequency = calculate_allele_freq(bp1_read_count, bp2_read_count, bp1_opposing_read_count, bp2_opposing_read_count, purity)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
