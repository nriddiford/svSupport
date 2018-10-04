
def classify_sv(bp1_sig, bp2_sig):

    for k1 in sorted(bp1_sig.iterkeys()):
        for k2 in sorted(bp2_sig.iterkeys()):
            if '_bp1' in k1:
                print "read before breakpoint 1"
                if 'bp2_' in k2:
                    print "read after breakpoint 2"
                    return 'DEL', '5to3'
                elif '_bp2' in k2:
                    print "read before breakpoint 2"
                    return 'BND', '3to3'
            elif 'bp1_' in k1:
                print "read after breakpoint 1"
                if 'bp2_' in k2:
                    print "read after breakpoint 2"
                    return 'BND', '5to5'
                elif '_bp2' in k2:
                    print "read before breakpoint 2"
                    return 'TANDUP', '3to5'

            print "%s: %s" % (k1, bp1_sig[k1])
            print "%s: %s" % (k2, bp2_sig[k2])


def classify_cnv(chrom, rdr):
    if chrom in ['X', 'Y']:
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
        elif (rdr > 1):
            cnv_type = 'Heterozygous duplication'

    print("%s on %s") % (cnv_type, chrom)
    return cnv_type