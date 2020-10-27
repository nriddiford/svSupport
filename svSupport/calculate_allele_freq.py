from __future__ import division
import math

class AlleleFrequency(object):
    def __init__(self, total_oppose, total_support, tumour_purity, chrom, sex):
        self.total_oppose = total_oppose
        self.total_support = total_support
        self.tumour_purity = tumour_purity
        self.chrom = chrom
        self.sex = sex


    def read_support_af(self):
        p = self.tumour_purity
        total_support = self.total_support
        total_oppose = self.total_oppose
        if total_support == 0: return 0, 0

        total_reads = total_support + total_oppose

        expected_oppose = (1-p)*total_reads
        # print("Expect to see %s/%s with purity at %s") % (expected_oppose, total_reads, p)
        adjusted_oppose = total_oppose-expected_oppose

        # print("Adjusted opposing reads: %s - %s = %s") % (total_oppose, expected_oppose, adjusted_oppose)

        if adjusted_oppose < 0:
            adjusted_oppose = 0
            print("Not sure if we should be here ...")

        allele_frequency = round(float(total_support)/(float(total_support)+float(total_oppose)), 2)

        if p == 1:
            adj_allele_frequency = allele_frequency
        else:
            adj_allele_frequency = float(total_support/( total_support + adjusted_oppose ))

        adj_allele_frequency = round(adj_allele_frequency, 2)
        print("* Allele frequency adjusted from %s to %s" % (allele_frequency, adj_allele_frequency))

        return allele_frequency, adj_allele_frequency


    def read_depth_af(self):
        n = self.total_oppose
        t = self.total_support
        p = self.tumour_purity
        c = self.chrom
        sex = self.sex

        su = abs(n - t)
        # print("Supporting reads = %s (%s-%s)") % (su, n, t)
        op = abs(n - su)
        # print("Opposing reads = %s (%s-%s)") % (op, n, su)

        r1 = round((t/n), 2)
        af = su/(op+su)

        adjop = op

        if p < 1:
            adjop = ((1-p) * op) + 0.01

        # print("Adjusting op reads for tp: %s = %s (%s*%s)") % (p, adjop, op, p)
        adjaf = su/(su+adjop)

        # print("%s / (%s + %s)") % (su, adjop, su)
        if r1 <= 1:
            r2 = round(t/(n + adjop), 2)
        else:
            r2 = round((t + adjop)/n, 2)

        if sex == 'XX':
            adjaf = round((adjaf/2), 2)
            af = round(af/2, 2)
        else:
            if c != 'X' and c != 'Y':
                adjaf = round((adjaf/2), 2)
                af = round(af/2, 2)
            else:
                adjaf = round(adjaf, 2)
                af = round(af, 2)

        log2_rd_ratio = round(math.log(r2, 2), 2)
        print("* Log2, purity-adjusted read depth ratio = %s " % log2_rd_ratio )
        print("* Allele frequency adjusted from %s to %s" % (af, adjaf))

        return af, adjaf, r2, log2_rd_ratio
