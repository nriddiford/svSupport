import os
import math
from itertools import islice
import pysam


def make_dirs(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def getChroms(options):
    """Read a file specifying native chromosomes"""
    chroms = {}
    with open(options.chromfile) as c:
        for line in c:
            chrom, length = line.strip().split()
            chroms[chrom] = length
            # chroms.append(str(line))
    return chroms


def cleanup(out_dir):
    print("Cleaning up old files in '%s'" % out_dir)
    out_dir = os.path.abspath(out_dir)
    for f in os.listdir(out_dir):
        extension = os.path.splitext(f)[1][1:]
        if extension in ['bam', 'bai']:
            try:
                abs_file = os.path.join(out_dir, f)
                os.remove(abs_file)
            except OSError:
                print("Can't remove %s" % abs_file)
                pass
    return(out_dir)


def print_options(bam_in, ratio, chrom, bp1, bp2, find_bps, debug, test, out_dir):
    options = ['Bam file', 'ratio', 'Chrom', 'bp1', 'bp2', 'hone_bps', 'debug', 'test', 'Out dir']
    args = [bam_in, ratio, chrom, bp1, bp2, find_bps, debug, test, out_dir]
    print("Running with options:")
    print("--------")
    for index, (value1, value2) in enumerate(zip(options, args)):
         print("o %s: %s") % (value1, value2)
    print("--------")
    print("python svSupport.py -i %s -l %s:%s-%s -f %s -t %s -d %d -o %s") % (bam_in, chrom, bp1, bp2, find_bps, test, debug, out_dir )
    print("--------")


def filterfn(read):
    """"Filter reads to ensure only properly paired, high quality reads are counted"""
    return (read.is_proper_pair and read.is_paired and read.tlen > 0 and not read.is_supplementary and not read.is_duplicate and not read.is_unmapped and not read.mate_is_unmapped)


def find_is_sd(bam_file, samplesize):
    """"Get empirical insert size distribution and return mean + 5 * SD"""
    bam = pysam.Samfile(bam_file, 'rb')
    l = islice((read.tlen for read in bam if filterfn(read)), samplesize)
    l = list(l)
    assert len(l) == samplesize
    mean = float(sum(l)) / len(l)
    sdev = math.sqrt(float(sum([(x - mean) ** 2 for x in l])) / (len(l) - 1))
    slop = int(mean + 5 * sdev)
    print('Using slop equal to 5 standard deviations from insert size mean: {:.0f}'.format(slop))
    return slop