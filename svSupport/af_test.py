def read_support_af(tumour_purity, total_support, total_oppose):
    total_reads = total_support + total_oppose

    expected_oppose = (1 - tumour_purity) * total_reads
    adjusted_oppose = total_oppose - expected_oppose

    if adjusted_oppose < 0:
        adjusted_oppose = 0

    adj_allele_frequency = float(total_support / (total_support + adjusted_oppose))
    adj_allele_frequency = round(adj_allele_frequency, 2)

    return adj_allele_frequency

read_support_af(0.67, 25, 11)