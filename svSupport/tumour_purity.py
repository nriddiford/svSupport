class Purity(object):
    def __init__(self, total_oppose, total_support, tumour_purity):

        self.total_oppose = total_oppose
        self.total_support = total_support
        self.tumour_purity = tumour_purity


    def get_af(self):
        p = self.tumour_purity
        total_support = self.total_support
        total_oppose = self.total_oppose

        total_reads = total_support + total_oppose
        expected_oppose = (1-p)*total_reads

        # print("Expect to see %s/%s with purity at %s") % (expected_oppose,total_reads, p)
        adjusted_oppose = total_oppose-expected_oppose

        # print("Adjusted opposing reads: %s - %s = %s") % (total_oppose, expected_oppose,adjusted_oppose)

        if adjusted_oppose < 0:
            adjusted_oppose = 0
            print("Not sure if we should be here ...")

        # print("Tumour purity set to %s" % p)
        allele_frequency = round(float(total_support)/(float(total_support)+float(total_oppose)),2)
        if p == 1:
            adj_allele_frequency = allele_frequency
        else:
            adj_allele_frequency = float( total_support/( total_support + adjusted_oppose ) )

        adj_allele_frequency = round(adj_allele_frequency, 2)
        print("Unadjusted allele frequency = %s" % allele_frequency)
        print("Adjusted allele frequency = %s" % adj_allele_frequency)
        print("------")


        return(adj_allele_frequency)
