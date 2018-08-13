import os
import math
from itertools import islice
import pysam


def classify_sv(bp1_best_guess, bp2_best_guess):
    if bp1_best_guess == 'F_bp1' and bp2_best_guess == 'bp2_R':
        sv_type = 'Deletion'
        read_sig = "__\ bp1 ... bp2 /__"
    elif bp1_best_guess == 'F_bp1' and bp2_best_guess == 'F_bp2':
        sv_type = 'Inversion type I'
        read_sig = "__\ bp1 ... __\ bp2"
    elif bp1_best_guess == 'bp1_R' and bp2_best_guess == 'bp2_R':
        sv_type = "Inversion type II"
        read_sig = "bp1 /__ ... bp2 /__"
    elif bp1_best_guess == 'bp1_R' and bp2_best_guess == 'F_bp2':
        sv_type = "Tandem duplication"
        read_sig = "bp1 /__ ... __\ bp2"
        bp1_best_guess, bp2_best_guess = 'F_bp1', 'bp2_R'
    else:
        sv_type = "Unknown - assuming deletion"
        read_sig = "__\ bp1 ... bp2 /__"
        bp1_best_guess, bp2_best_guess = 'F_bp1', 'bp2_R'

    return(bp1_best_guess, bp2_best_guess, sv_type, read_sig)


def classify_cnv(chrom, rdr):
    if chrom == 'X' or chrom == 'Y':
        if (rdr < 1):
            cnv_type = 'Homozygous deletion'
        elif (rdr >= 3):
            cnv_type = 'Homozygous triplication'
        elif (rdr >= 1):
            cnv_type = 'Homozygous duplication'
    else:
        if (rdr <= 0.5):
            cnv_type = 'Homozygous deletion'
        elif (rdr <= 1):
            cnv_type = 'Heterozygous deletion'
        elif (rdr >= 2):
            cnv_type = 'Heterozygous triplication'
        elif (rdr >= 1.5):
            cnv_type = 'Heterozygous duplication'

    print("%s on %s") % (cnv_type, chrom)
    return cnv_type


def make_dirs(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def getChroms(options):
    """Read a file specifying native chromosomes"""
    chroms = []
    with open(options.chromfile) as c:
        for line in c:
            line = line.strip()
            chroms.append(str(line))
    return chroms


def cleanup(out_dir):
    print("Cleaning up old files in '%s'" % out_dir)
    out_dir = os.path.abspath(out_dir)
    for f in os.listdir(out_dir):
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