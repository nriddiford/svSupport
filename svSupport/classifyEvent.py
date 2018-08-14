
def classify_sv1(bp1_best_guess, bp2_best_guess):
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



def classify_sv(bp_sig):



    return bp_sig


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