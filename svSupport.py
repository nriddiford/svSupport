import sys, os, re
import pysam
from optparse import OptionParser

out_dir = 'out'

def bp1_supporting_reads(bamFile, chrom, bp1, bp2, slop):
    samfile = pysam.Samfile(bamFile, "rb")
    start=bp1-slop
    bp1_reads = []
    out_file = os.path.join(out_dir, "bp1_sv_reads" + ".bam")
    bp1_sv_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
    count = 0
    for read in samfile.fetch(chrom, start, bp1):
        read_end_pos = read.pos + read.alen
        mate_end_pos = read.mpos + read.alen

        if not read.is_proper_pair and not read.is_reverse and read.mpos > bp2:
            print("* bp1 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.pos, read_end_pos, read.mpos, mate_end_pos)
            bp1_sv_reads.write(read)
            bp1_reads.append(read.qname)
            count += 1

        if read_end_pos == bp1:
            print("* bp1 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.pos, read_end_pos)
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
        read_end_pos = read.pos + read.alen
        mate_end_pos = read.mpos + read.alen

        if not read.is_proper_pair and read.is_reverse and mate_end_pos < bp1:
            print("* bp2 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.pos, read_end_pos, read.mpos, mate_end_pos)
            bp2_sv_reads.write(read)
            bp2_reads.append(read.qname)
            count += 1

        # mapped_start = read_end_pos - read.alen

        if read.pos +1 == bp2:
            print("* bp2 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.pos, read_end_pos)
            bp2_reads.append(read.qname)
            bp2_sv_reads.write(read)
            count += 1

    bp2_sv_reads.close()
    pysam.index(out_file)

    return(bp2_reads, count)


def merge_bams(out_file, bams):
    out_file = os.path.join(out_dir, out_file)
    in_files = ', '.join(bams)
    print("Merging bam files '%s' into '%s'") % (in_files, out_file)
    merge_parameters = ['-f', out_file] + bams
    pysam.merge(*merge_parameters)
    pysam.index(out_file)


def findSupport(bam_in, chrom, bp1, bp2, slop):
    bp1_sv_reads, bp1_read_count = bp1_supporting_reads(bam_in, chrom, bp1, bp2, slop)
    bp2_sv_reads, bp2_read_count = bp2_supporting_reads(bam_in, chrom, bp1, bp2, slop)
    print(bp2_sv_reads)

    to_merge = ["out/bp1_sv_reads.bam", "out/bp2_sv_reads.bam"]
    merge_bams("sv_support.bam", to_merge)
    print(bp1_read_count, bp2_read_count)
    return(bp1_sv_reads, bp1_read_count, bp2_sv_reads, bp2_read_count)


def main():
    parser = OptionParser()

    parser.add_option("-i", \
                      "--inFile", \
                      dest="inFile",
                      action="store",
                      help="A sorted .bam file containing the reads " + \
                           "supporting the structural variant calls", \
                           metavar="FILE")

    parser.add_option("-s", \
                      "--slop", \
                      dest="slop",
                      type="int",
                      default=500,
                      action="store",
                      help="Distance from breakpoint to look for reads" + \
                           "Default: 500")

    parser.add_option("-l", \
                      "--loci", \
                      dest="region",
                      type="string",
                      action="store",
                      help="The chromosome and breakpoints for a " + \
                           "structural variant in the format: " + \
                           "'chrom:bp_1-bp_2'")

    options, args = parser.parse_args()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if options.inFile is None:
        parser.print_help()
        print
    else:
        try:
            bam_in = options.inFile
            region = options.region
            slop   = options.slop

            chrom, bp1, bp2 = re.split(':|-', region)
            bp1 = int(bp1)
            bp2 = int(bp2)
            print("-----\nBam file: '%s'\nChrom: %s\nbp1: %s\nbp2: %s\nslop: %s\n-----") % (bam_in, chrom, bp1, bp2, slop)
            findSupport(bam_in, chrom, bp1, bp2, slop)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
