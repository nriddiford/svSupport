from optparse import OptionParser


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

    parser.add_option("--chromosomes",
                      dest="chromfile",
                      action="store",
                      help="A file listing chromosome names to consider for normal mapping",
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

    parser.add_option("-s",
                      "--slop",
                      dest="slop",
                      action="store",
                      type="int",
                      help="Explicitly set the distance from breakpoints " +
                           "to consider as informative for SV")

    parser.add_option("-l",
                      "--loci",
                      dest="region",
                      action="store",
                      help="The chromosome and breakpoints for a " +
                           "structural variant in the format: " +
                           "'chrom:bp_1-bp_2' or 'chrom1:bp_1-chrom2:bp_2")
    parser.add_option("--sex",
                      dest="sex",
                      action="store",
                      help="Sex of individual: XX, XY [Default: XY] ")

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

    parser.set_defaults(out_dir='out',
                        purity=1,
                        chromfile='chrom_lengths.txt')

    options, args = parser.parse_args()

    if (options.in_file is None or options.region is None) and not options.test and not options.config:
        parser.print_help()
        print

    return (options, args)