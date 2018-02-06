#!/usr/bin/env python
from __future__ import division
import os, re
import pysam
import sys
sys.dont_write_bytecode = True

from optparse import OptionParser
from deletions import Deletions


def calculate_allele_freq(bp1_read_count, bp2_read_count, bp1_opposing_read_count, bp2_opposing_read_count, tumour_purity):
    total_support =  bp1_read_count + bp2_read_count
    total_oppose = bp1_opposing_read_count + bp2_opposing_read_count

    print("Tumour purity set to %s" % tumour_purity)

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


def print_options(bam_in, chrom, bp1, bp2, slop, find_bps, debug, test, out_dir):
    options = ['Bam file', 'Chrom', 'bp1', 'bp2', 'slop', 'search_bps', 'debug', 'test', 'Out dir']
    args = [bam_in, chrom, bp1, bp2, slop, find_bps, debug, test, out_dir]
    print("Running with options:")
    print("--------")
    for index, (value1, value2) in enumerate(zip(options, args)):
         print("o %s: %s") % (value1, value2)
    print("--------")
    print("python svSupport.py -i %s -l %s:%s-%s -s %s -f %s -t %s -d %d -o %s") % (bam_in, chrom, bp1, bp2, slop, find_bps, test, debug, out_dir )
    print("--------")

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

    parser.add_option("-t", \
                    "--test", \
                    dest="test",
                    action="store_true",
                    help="Run on test data")

    parser.set_defaults(slop=500, out_dir='../out', debug=0, test=False, purity=1, find_bps=0)
    options, args = parser.parse_args()

    if (options.in_file is None or options.region is None) and options.test is False:
      parser.print_help()
      print

    return(options, args)


def main():
    options, args = get_args()

    global out_dir
    global debug

    if options.test is True:
        print
        print("Running in test mode...")
        print

        options.region = '3L:9892365-9894889'
        options.out_dir = '../test_out'
        options.in_file = '../data/test.bam'
        options.debug = True

    if options.in_file is not None and options.region is not None:
        try:
            bam_in = options.in_file
            region = options.region
            slop   = options.slop
            out_dir = options.out_dir
            debug = options.debug
            purity = float(options.purity)
            find_bps = options.find_bps
            test = options.test

            chrom, bp1, bp2 = re.split(':|-', region)
            bp1 = int(bp1)
            bp2 = int(bp2)
            slop = int(slop)

            if debug:
                print_options(bam_in, chrom, bp1, bp2, slop, find_bps, debug, test, out_dir)

            #-------------------
            # Adjust breakpoints
            #-------------------

            if find_bps:
                bp1, bp1_count = search_bps(bam_in, chrom, bp1, 'bp1')
                bp2, bp2_count = search_bps(bam_in, chrom, bp2, 'bp2')
                print("Bp1 adjusted to: %s [%s split reads found]") % (bp1, bp1_count)
                print("Bp2 adjusted to: %s [%s split reads found]") % (bp2, bp2_count)

            make_dirs(bam_in, out_dir)

            #-------------------
            # DELETIONS
            #-------------------

            del_support = Deletions(bam_in, chrom, bp1, bp2, slop, 'support', out_dir, debug)
            del_oppose  = Deletions(bam_in, chrom, bp1, bp2, slop, 'oppose', out_dir, debug)

            bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count = del_support.get_reads()
            bp1_opposing_reads, bp1_opposing_read_count, bp2_opposing_reads, bp2_opposing_read_count = del_oppose.get_reads()

            #-------------------
            # Calculate af
            #-------------------

            allele_frequency = calculate_allele_freq(bp1_read_count, bp2_read_count, bp1_opposing_read_count, bp2_opposing_read_count, purity)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
