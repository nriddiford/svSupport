import operator


def classify_sv(bp1_sig, bp2_sig):
    most_common_bp1 = sorted(bp1_sig.items(), key=operator.itemgetter(1), reverse=True)
    most_common_bp2 = sorted(bp2_sig.items(), key=operator.itemgetter(1), reverse=True)

    for k1, v1 in most_common_bp1:
        for k2, v2 in most_common_bp2:
            if '_bp1' in k1:
                print("read before breakpoint 1")
                if 'bp2_' in k2:
                    print("read after breakpoint 2")
                    return 'DEL', '5to3', k1, k2
                elif '_bp2' in k2:
                    print("read before breakpoint 2")
                    return 'BND', '3to3', k1, k2
            elif 'bp1_' in k1:
                print("read after breakpoint 1")
                if 'bp2_' in k2:
                    print("read after breakpoint 2")
                    return 'BND', '5to5', k1, k2
                elif '_bp2' in k2:
                    print("read before breakpoint 2")
                    return 'TANDUP', '3to5', k1, k2


def classify_cnv(chrom, rdr, sex):
    if chrom in ['X', 'Y'] and sex == 'XY':
        if (rdr < 1):
            cnv_type = 'Homozygous deletion'
        elif (rdr >= 3):
            cnv_type = 'Homozygous triplication'
        elif (rdr >= 1):
            cnv_type = 'Homozygous duplication'
    else:
        if (rdr <= 0.3):
            cnv_type = 'Homozygous deletion'
        elif (rdr <= 1):
            cnv_type = 'Heterozygous deletion'
        elif (rdr >= 2):
            cnv_type = 'Heterozygous triplication'
        elif (rdr > 1):
            cnv_type = 'Heterozygous duplication'

    print("%s on %s" % (cnv_type, chrom))
    return cnv_type