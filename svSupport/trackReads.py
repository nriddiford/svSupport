"""

#------------------------
# pySam variables
#------------------------

o next_reference_start   =   the position of the mate/next read (formerly 'mpos')
o reference_start        =   0-based leftmost coordinate (formerly 'pos')
o reference_length       =   aligned length of the read on the reference genome
                             This is equal to aend - pos. Returns None if not available (formerly 'alen')
o query_name             =   Read name (formaerly 'qname')

"""

class TrackReads(object):
    def __init__(self, read, mate, chrom1, chrom2, duplicates):

        self.read = read
        self.mate = mate
        self.chrom1 = chrom1
        self.chrom2 = chrom2
        self.duplicates = duplicates


    def check_for_standard_dup(self):
        dupkey = '_'.join(map(str, [self.read.reference_start, self.read.reference_end, self.mate.next_reference_start]))
        dup = False
        self.duplicates[dupkey] += 1
        if self.duplicates[dupkey] > 1:
            dup = True

        return self.duplicates, dup


    def check_for_disc_dup(self):
        dupkey = '_'.join(map(str, [self.chrom1, self.read.reference_start, self.read.reference_end, self.chrom2, self.mate.next_reference_start]))
        dup = False
        self.duplicates[dupkey] += 1
        if self.duplicates[dupkey] > 1:
            dup = True

        return self.duplicates, dup


    def check_for_clipped_dup(self):
        dup = False
        try:
            self.read.get_tag('SA')
            dupkey = '_'.join(map(str, [self.read.reference_start, self.read.reference_end, self.read.get_tag('SA').split(',')[1]]), self.chrom2, self.mate.next_reference_start)
            self.duplicates[dupkey] += 1
            if self.duplicates[dupkey] > 1:
                dup = True
        except:
            pass

        return self.duplicates, dup