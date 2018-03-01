# from svSupport.svSupport import classify_cnv

class Allele_frequency(object):
    def __init__(self, total_oppose, total_support, tumour_purity, read_depth_ratio, chrom):
        self.total_oppose = total_oppose
        self.total_support = total_support
        self.tumour_purity = tumour_purity
        self.read_depth_ratio = read_depth_ratio
        self.chrom = chrom


    def read_support_af(self):
        p = self.tumour_purity
        total_support = self.total_support
        total_oppose = self.total_oppose

        total_reads = total_support + total_oppose

        expected_oppose = (1-p)*total_reads
        # print("Expect to see %s/%s with purity at %s") % (expected_oppose, total_reads, p)
        adjusted_oppose = total_oppose-expected_oppose

        # print("Adjusted opposing reads: %s - %s = %s") % (total_oppose, expected_oppose, adjusted_oppose)

        if adjusted_oppose < 0:
            adjusted_oppose = 0
            print("Not sure if we should be here ...")

        # print("Tumour purity set to %s" % p)

        allele_frequency = round(float(total_support)/(float(total_support)+float(total_oppose)),2)
        if p == 1:
            adj_allele_frequency = allele_frequency
        else:
            adj_allele_frequency = float(total_support/( total_support + adjusted_oppose ))

        adj_allele_frequency = round(adj_allele_frequency, 2)


        print("Unadjusted allele frequency = %s" % allele_frequency)
        print("Adjusted allele frequency = %s" % adj_allele_frequency)
        return(adj_allele_frequency)


    def read_depth_af(self):
        n = self.total_oppose
        t = self.total_support
        p = self.tumour_purity
        c = self.chrom

        su = abs(n - t)
        # print("Supporting reads = %s (%s-%s)") % (su, n, t)
        op = abs(n - su)
        r1 = round((t/n), 2)
        print("* read depth ratio = %s " % r1)
        af = su/(op+su)

        adjop = (op ** p) + 0.01
        adjaf = su/(su+adjop)
        # print("%s / (%s + %s)") % (su, adjop, su)
        if r1 <=1:
            r2 = round(t/(n+(n*p)), 2)
            print("* Purity-adjusted read depth ratio = %s " % r2)
        else:
            r2 = round((t+(t*p))/n, 2)
            print("* Purity-adjusted read depth ratio = %s " % r2)

        if c != 'X' and c != 'Y':
            adjaf = round((adjaf/2), 2)
            af = round(af/2, 2)
        else:
            adjaf = round(adjaf, 2)
            af = round(af, 2)

        print("* allele frequency = %s " % af)
        print("* Purity-adjusted allele frequency = %s " % adjaf)
        return(adjaf, r2)