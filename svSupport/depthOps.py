from __future__ import division
import pysam

def get_depth(bam_in, normal, chrom, bp1, bp2, chroms):
    """Get the number of mapped reads in both t and n bams accross all chroms
       Then get the number of mapped reads within the CNV region
       """
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = [chrom]

    t_reads_by_chrom, tumour_mapped = count_reads(bam_in, chromosomes)
    t_read_count = region_depth(bam_in, chrom, bp1, bp2)

    n_reads_by_chrom, normal_mapped = count_reads(normal, chromosomes)
    n_read_count = region_depth(normal, chrom, bp1, bp2)

    """Mark as FP CN events where #reads/length < fraction"""
    # length = bp2 - bp1
    # v = n_read_count/length

    mapped_ratio = tumour_mapped / normal_mapped

    # Changed from mapped_ratio < 1 - 3.8.18
    if mapped_ratio >= 1:
        t_corr = t_read_count
        n_corr = int(round((n_read_count * mapped_ratio)))
    else:
        n_corr = n_read_count
        t_corr = round((t_read_count * mapped_ratio))

    adj_ratio = round((t_corr / n_corr), 2)

    return (n_corr, t_corr, adj_ratio)


def count_reads(bamfile, chromosomes):
    """Count the total number of mapped reads in a BAM file, filtering
    for chromosomes in `chromosomes`
    """

    lines = pysam.idxstats(bamfile)

    if type(lines) is str:
        lines = lines.strip().split('\n')

    total_mapped = 0
    total_mapped_chrom = {}

    for line in lines:
        chrom, len, mapped, unmapped = line.split('\t')
        if chrom in chromosomes:
            total_mapped_chrom[chrom] = int(mapped)
            total_mapped += int(mapped)

    print("Total number of mapped reads on chroms %s %s: %s") % (chromosomes, bamfile, total_mapped)
    return(total_mapped_chrom, total_mapped)


def region_depth(bamfile, chrom, bp1, bp2):
    """Count the total number of mapped reads in a genomic region"""

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
