import pysam
from collections import defaultdict
import math

def get_depth(bam_in, normal, chrom, bp1, bp2):
    """Get the number of mapped reads in both t and n bams accross all chroms
       Then get the number of mapped reads within the CNV region

       """
    chromosomes = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
    t_reads_by_chrom, tumour_mapped = count_reads(bam_in, chromosomes)
    t_read_count = region_depth(bam_in, chrom, bp1, bp2)


    n_reads_by_chrom, normal_mapped = count_reads(normal, chromosomes)
    n_read_count = region_depth(normal, chrom, bp1, bp2)

    mapped_ratio = tumour_mapped / normal_mapped

    print ("t: %s" % tumour_mapped)
    print ("n: %s" % normal_mapped)

    # Changed from mapped_ratio < 1 - 3.8.18
    if mapped_ratio >= 1:
        t_corr = t_read_count
        n_corr = round((n_read_count * mapped_ratio))
    else:
        n_corr = n_read_count
        t_corr = round((t_read_count * mapped_ratio))

    adj_ratio = round((t_corr / n_corr), 2)

    print math.log(adj_ratio, 2)

    print("Normalised read count ratio: %s (%s/%s)") % (adj_ratio, t_corr, n_corr)
    return (n_corr, t_corr, adj_ratio)


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
