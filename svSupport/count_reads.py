import pysam
from collections import defaultdict

def count_reads(bamfile, chromosomes):
    """Count the total number of mapped reads in a BAM file, filtering
    the chromosome given in chroms_to_include list
    """
    lines = pysam.idxstats(bamfile)

    if type(lines) is str:
        lines = lines.strip().split('\n')

    total_mapped = 0
    total_mapped_chrom = defaultdict(int)

    for line in lines:
        chrom, _len, nmapped, _nunmapped = line.split('\t')
        if chrom in chromosomes:
            total_mapped_chrom[chrom] += int(nmapped)
            total_mapped += int(nmapped)

    print("Total number of mapped reads in %s: %s") % (bamfile, total_mapped)
    return(total_mapped_chrom, total_mapped)

def region_depth(bamfile, chrom, bp1, bp2):
    """Count the total number of mapped reads in a genomic region, filtering
    the chromosome given in chroms_to_ignore list
    """

    samfile = pysam.Samfile(bamfile, "rb")
    count = 0
    for read in samfile.fetch(chrom, bp1, bp2):
        if read.is_unmapped:
            continue
        if read.mapq < 3:
            continue

        count += 1

    print("Reads in %s:%s-%s: %s") % (chrom, bp1, bp2, count)
    return(count)

def scale_factor(mappped_reads):
    scale_factor = 1
    # the CPM (norm is based on post-filtering total counts of reads in BAM "bam_mapped")
    million_reads_mapped = float(mappped_reads) / 1e6
    scale_factor *= 1.0 / (million_reads_mapped)
    print("scale_factor = %s" % scale_factor)
