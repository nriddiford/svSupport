class Purity(object):
    def __init__(self, total_oppose, total_support, tumour_purity, read_depth_ratio):

        self.total_oppose = total_oppose
        self.total_support = total_support
        self.tumour_purity = tumour_purity
        self.read_depth_ratio = read_depth_ratio

    def get_af(self):
        p = self.tumour_purity
        total_support = self.total_support
        total_oppose = self.total_oppose

        total_reads = total_support + total_oppose
        expected_oppose = (1-p)*total_reads

        if self.read_depth_ratio:
            if self.read_depth_ratio < 1:
                print ("Deletion")
                adjusted_oppose = expected_oppose+(expected_oppose*p)
                adjusted_ratio = total_support/adjusted_oppose
                print(adjusted_ratio)
                allele_frequency = 1 - adjusted_ratio
            elif self.read_depth_ratio <= 2:
                print ("Duplication")
                expected_support = (1-p)*total_reads
                adjusted_support = expected_support+(expected_support*p)
                adjusted_ratio = adjusted_support/total_oppose
                print(adjusted_ratio)
                allele_frequency = 1 - (2 - adjusted_ratio)

            elif self.read_depth_ratio <= 3:
                print ("Triplication")
                expected_support = (1-p)*total_reads
                adjusted_support = expected_support+(expected_support*p)
                adjusted_ratio = adjusted_support/total_oppose
                print(adjusted_ratio)
                allele_frequency = 1 - (3 - adjusted_ratio)

            return(round(allele_frequency,2))

        print("Total oppose in region: %s " % total_oppose)
        print("Total support in region: %s " % total_support)

        print("Expect to see %s/%s with purity at %s") % (expected_oppose, total_reads, p)
        adjusted_oppose = total_oppose-expected_oppose

        print("Adjusted opposing reads: %s - %s = %s") % (total_oppose, expected_oppose, adjusted_oppose)

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
