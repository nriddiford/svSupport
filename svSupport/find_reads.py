import pysam
import os, re, sys
from merge_bams import merge_bams

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

class FindReads(object):
    def __init__(self, bam_in, chrom, bp1, bp2, slop, out_dir, debug, bp1_class, bp2_class):


        self.bam_in = bam_in
        self.chrom = chrom
        self.bp1 = bp1
        self.bp2 = bp2
        self.slop = slop
        self.out_dir = out_dir
        self.debug = debug
        self.bp1_class = bp1_class
        self.bp2_class = bp2_class

    def set_window(self, breakpoint):

        if self.bp1_class == 'F_bp1':
            bp1_start = self.bp1 - self.slop
            bp1_end = self.bp1
        elif self.bp1_class == 'bp1_R':
            bp1_start = self.bp1
            bp1_end = self.bp1 + self.slop

        if self.bp2_class == 'bp2_R':
            bp2_start = self.bp2
            bp2_end = self.bp2 + self.slop
        elif self.bp2_class == 'F_bp2':
            bp2_start = self.bp2 - self.slop
            bp2_end = self.bp2

        if breakpoint == 'bp1':
            return(bp1_start, bp1_end)
        else:
            return(bp2_start, bp2_end)


    def bp1_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        bp1_start, bp1_end = self.set_window('bp1')
        supporting_reads = []
        opposing_reads = []

        support_out = os.path.join(self.out_dir, "bp1_supporting_reads" + ".bam")
        oppose_out = os.path.join(self.out_dir, "bp1_opposing_reads" + ".bam")

        with pysam.AlignmentFile(support_out, "wb", template=samfile) as bp1_supporting_reads, pysam.AlignmentFile(oppose_out, "wb", template=samfile) as bp1_opposing_reads:
            for read in samfile.fetch(self.chrom, bp1_start, bp1_end):
                read_end_pos = read.reference_start + read.reference_length
                try:
                    mate = samfile.mate(read)
                    mate_end_pos = mate.reference_start + mate.reference_length
                except:
                     continue

                if read.is_duplicate:
                    continue

                read_start = read.reference_start +1
                mate_start = mate.reference_start +1

                if read.reference_name != mate.reference_name:
                    continue

                # deletion
                if self.bp1_class == 'F_bp1' and self.bp2_class == 'bp2_R':
                    # supporting
                    if read_end_pos <= self.bp1 and mate_start >= self.bp2:
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)
                    elif read_end_pos == self.bp1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if (read_end_pos < self.bp1 and mate_start > self.bp1 and mate_start < self.bp2) or (read.reference_start > self.bp1 and mate_end_pos < self.bp1):
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        elif read.reference_start < self.bp1 and read_end_pos > self.bp1:
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                # type I inversion - ('F_bp1', 'F_bp2')
                elif self.bp1_class == 'F_bp1' and self.bp2_class == 'F_bp2':
                    if read.reference_name != mate.reference_name:
                        continue
                    # supporting
                    if not read.is_proper_pair and read_end_pos <= self.bp1 and mate_end_pos >= self.bp1 and mate_end_pos <= self.bp2:
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)

                    elif read_end_pos == self.bp1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if (read_end_pos < self.bp1 and mate_start > self.bp1 and mate_start < self.bp2) or (read.reference_start > self.bp1 and mate_end_pos < self.bp1):
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        elif read.reference_start < self.bp1 and read_end_pos > self.bp1:
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                # type II inversion - ('bp1_R', 'bp2_R')
                elif self.bp1_class == 'bp1_R' and self.bp2_class == 'bp2_R':
                    if read.reference_name != mate.reference_name:
                        continue
                    # supporting
                    if not read.is_proper_pair and read_start >= self.bp1 and read_start <= self.bp2 and mate_start >= self.bp2:
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)

                    # change 13.2.18 - need to sort out 0/1 based start - this works for type II...
                    elif read_start == self.bp1 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp1(read, mate, bp1_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if ( read_start > self.bp1 and mate_start < self.bp1 ) or (read_end_pos < self.bp1 and mate_start > self.bp1):
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        elif read_start < self.bp1 and read_end_pos > self.bp1:
                            self.print_and_write_bp1(read, mate, bp1_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                else:
                    print("Reads must be classified")
                    sys.exit()

        support_count = len(set(supporting_reads))
        oppose_count = len(set(opposing_reads))

        return(supporting_reads, support_count, support_out, opposing_reads, oppose_count, oppose_out)


    def print_and_write_bp1(self, read, mate, out_bam, read_list, evidence, for_against, read_end_pos, mate_end_pos):
        if self.debug:
            print("* bp1 %s %s read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (evidence, for_against, read.query_name, read.seq, read.reference_start, read_end_pos, mate.reference_start, mate_end_pos)
        if read.query_name not in read_list:
            out_bam.write(read)
            if evidence == 'mate_pair' or evidence == 'disc_read':
                out_bam.write(mate)
            read_list.append(read.query_name)


    def bp2_reads(self, bp1_supporting_reads, bp1_opposing_reads):
        samfile = pysam.Samfile(self.bam_in, "rb")
        bp2_start, bp2_end = self.set_window('bp2')
        supporting_reads = []
        opposing_reads = []

        support_out = os.path.join(self.out_dir, "bp2_supporting_reads" + ".bam")
        oppose_out = os.path.join(self.out_dir, "bp2_opposing_reads" + ".bam")

        with pysam.AlignmentFile(support_out, "wb", template=samfile) as bp2_supporting_reads, pysam.AlignmentFile(oppose_out, "wb", template=samfile) as bp2_opposing_reads:
            for read in samfile.fetch(self.chrom, bp2_start, bp2_end):
                read_end_pos = read.reference_start + read.reference_length
                try:
                    mate = samfile.mate(read)
                    mate_end_pos = mate.reference_start + mate.reference_length
                except:
                     continue

                if read.is_duplicate:
                    continue

                read_start = read.reference_start +1
                mate_start = mate.reference_start +1

                if read.reference_name != mate.reference_name:
                    continue

                # deletion
                if self.bp1_class == 'F_bp1' and self.bp2_class == 'bp2_R':
                    # supporting
                    if read_start >= self.bp2 and mate_end_pos <= self.bp1:
                    # if not read.is_proper_pair and read.is_reverse and mate_end_pos < self.bp1:
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)
                    elif read_start == self.bp2 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if (read_end_pos < self.bp2 and mate_start > self.bp2) or (read_start > self.bp2 and mate_end_pos < self.bp2 and mate_end_pos > self.bp1):
                        # if (read_end_pos < self.bp1 and mate.reference_start > self.bp1 and mate.reference_start < self.bp2) or (read.reference_start > self.bp1 and mate_end_pos < self.bp1):
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        # elif read.reference_start < self.bp1 and read_end_pos > self.bp1:
                        elif read_start < self.bp2 and read_end_pos > self.bp2:
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                # type I inversion - ('F_bp1', 'F_bp2')
                elif self.bp1_class == 'F_bp1' and self.bp2_class == 'F_bp2':
                    # supporting
                    if not read.is_proper_pair and read_end_pos +1 <= self.bp2 and mate_end_pos <= self.bp1:
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)

                    elif read_end_pos == self.bp2 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if (read_end_pos < self.bp2 and mate_start > self.bp2) or (read_start > self.bp2 and mate_end_pos < self.bp2):
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        elif read_start < self.bp2 and read_end_pos > self.bp2:
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                # type II inversion - ('bp1_R', 'bp2_R')
                elif self.bp1_class == 'bp1_R' and self.bp2_class == 'bp2_R':
                    # supporting
                    if not read.is_proper_pair and read_start >= self.bp2 and mate_start >= self.bp1 and mate_end_pos < self.bp2:
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'disc_read', 'supporting', read_end_pos, mate_end_pos)

                    elif read_start == self.bp2 and re.findall(r'(\d+)[S|H]', read.cigarstring):
                        self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_supporting_reads, supporting_reads, 'clipped_read', 'supporting', read_end_pos, mate_end_pos)

                    # opposing
                    elif read.qname not in supporting_reads:
                        if ( read_end_pos < self.bp2 and mate_start > self.bp2 ) or (read_start > self.bp2 and mate_end_pos < self.bp2):
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'mate_pair', 'opposing', read_end_pos, mate_end_pos)
                        elif read_start < self.bp2 and read_end_pos > self.bp2:
                            self.print_and_write_bp2(read, mate, bp1_supporting_reads, bp1_opposing_reads, bp2_opposing_reads, opposing_reads, 'spanning', 'opposing', read_end_pos, mate_end_pos)

                else:
                    print("Reads must be classified")
                    sys.exit()
        support_count = len(set(supporting_reads))
        oppose_count = len(set(opposing_reads))

        return(supporting_reads, support_count, support_out, opposing_reads, oppose_count, oppose_out)


    def print_and_write_bp2(self, read, mate, bp1_supporting_reads, bp1_opposing_reads, out_bam, read_list, evidence, for_against, read_end_pos, mate_end_pos):
        if self.debug:
            print("* bp2 %s %s read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (evidence, for_against, read.query_name, read.seq, read.reference_start +1, read_end_pos, mate.reference_start +1, mate_end_pos)
        if read.query_name not in bp1_supporting_reads and read.query_name not in bp1_opposing_reads and read.query_name not in read_list:
            out_bam.write(read)
            if evidence == 'mate_pair' or evidence == 'disc_read':
                out_bam.write(mate)
            read_list.append(read.query_name)
