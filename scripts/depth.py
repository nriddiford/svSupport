import pysam
from collections import defaultdict

filename='../data/A373R13.tagged.filt.SC.RG.bam'

def bedtools(filename):
    """simulate the behaviour of bedtools"""
    bamfile = pysam.AlignmentFile(filename, 'rb')
    read_count = defaultdict(int)

    for ref in bamfile.header['SQ']:
        name = ref['SN']
        pileup = bamfile.pileup()
        for pos, column in enumerate(pileup, 1):
            depth = column.nsegments
            read_count[name] += depth

        for c in read_count:
            print(c, read_count[c])

# bedtools(filename)

chroms_to_include = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
def count_reads(filename, chroms_to_include):
    """Count the total number of mapped reads in a BAM file, filtering
    the chromosome given in chroms_to_ignore list
    """
    lines = pysam.idxstats(filename)

    if type(lines) is str:
        lines = lines.strip().split('\n')

    total_mapped = 0
    total_mapped_chrom = defaultdict(int)

    for line in lines:
        chrom, _len, nmapped, _nunmapped = line.split('\t')
        if chrom in chroms_to_include:
            total_mapped_chrom[chrom] += int(nmapped)
            total_mapped += int(nmapped)

    for c in total_mapped_chrom:
        print(c, total_mapped_chrom[c])

    return(total_mapped_chrom, total_mapped)

count_reads(filename, chroms_to_include)


def region_depth(filename, chrom, bp1, bp2):
    """Count the total number of mapped reads in a genomic region, filtering
    the chromosome given in chroms_to_ignore list
    """

    samfile = pysam.Samfile(filename, "rb")
    count = 0
    for read in samfile.fetch(chrom, bp1, bp2):
        count += 1

    return(count)
region_depth(filename, 'X', 3059398,3145934)
