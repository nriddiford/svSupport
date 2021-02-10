from __future__ import division
import pysam
from getReads import filterContamination


def get_depth(bam_in, normal, chrom, bp1, bp2, chroms, notes, options, chrom_dict):
    """Get the number of mapped reads in both t and n bams accross all chroms
       Then get the number of mapped reads within the CNV region
       """
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = [chrom]

    t_reads_by_chrom, tumour_mapped = count_reads(bam_in, chromosomes)
    t_read_count, t_contamination_count, read_length = region_depth(bam_in, chrom, bp1, bp2, options)

    n_reads_by_chrom, normal_mapped = count_reads(normal, chromosomes)
    n_read_count, n_contamination_count, read_length = region_depth(normal, chrom, bp1, bp2, options)

    av_depth = n_reads_by_chrom[chrom]*read_length/int(chrom_dict[chrom])

    length = bp2 - bp1
    v = (n_read_count/length)*read_length
    # print("Av depth in region", v, av_depth)

    if v/av_depth < 0.5 and v < 15:
        print(v, av_depth, v/av_depth, read_length)
        notes.append("low depth in normal bam")

    mapped_ratio = tumour_mapped / normal_mapped

    # Changed from mapped_ratio < 1 - 3.8.18
    # if mapped_ratio >= 1:
    #     t_corr = t_read_count
    #     n_corr = int(round((n_read_count * mapped_ratio)))
    # else:
    #     n_corr = n_read_count
    #     t_corr = int(round((t_read_count * mapped_ratio)))

    t_corr = t_read_count
    n_corr = int(round((n_read_count * mapped_ratio)))

    adj_ratio = round((t_corr / n_corr), 2)

    t_contamination_fraction = (t_contamination_count + 0.01)/t_read_count
    n_contamination_fraction = (n_contamination_count + 0.01)/n_read_count

    contamination_ratio = (t_contamination_fraction + 0.01) / (n_contamination_fraction + 0.01)

    if t_contamination_fraction >= 0.1 and contamination_ratio >= 2 and adj_ratio > 1:
        note = ' '.join(map(str, ["t_contamination:", t_contamination_count, "n_contamination:", n_contamination_count]))
        notes.append(note)

    return n_corr, t_corr, adj_ratio, notes


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
        chrom, length, mapped, unmapped = line.split('\t')
        if chrom in chromosomes:
            total_mapped_chrom[chrom] = int(mapped)
            total_mapped += int(mapped)

    print("Total number of mapped reads on chroms %s %s: %s" % (chromosomes, bamfile, total_mapped))
    return total_mapped_chrom, total_mapped


def region_depth(bamfile, chrom, bp1, bp2, options):
    """Count the total number of mapped reads in a genomic region"""

    samfile = pysam.Samfile(bamfile, "rb")
    count = 0
    contamination_count = 0
    read_lengths = 0

    check_read_length = 0

    for read in samfile.fetch(chrom, bp1, bp2):
        if read.is_unmapped:
            continue
        if read.mapq < 3:
            continue

        contaminated_read, contaminated_at_bp = filterContamination(read, bp1, options)
        if contaminated_read:
            contamination_count += 1
            continue

        if check_read_length < 100:
            read_lengths += read.infer_read_length()
            check_read_length += 1

        count += 1

    av_read_length = read_lengths/check_read_length
    print("Reads in %s:%s-%s: %s" % (chrom, bp1, bp2, count))
    return count, contamination_count, av_read_length
