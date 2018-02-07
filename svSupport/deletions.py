import pysam
import os, re
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

class Deletions(object):
    def __init__(self, bam_in, chrom, bp1, bp2, slop, read_type, out_dir, supporting_reads, debug=False):

        self.bam_in = bam_in
        self.chrom = chrom
        self.bp1 = bp1
        self.bp2 = bp2
        self.slop = slop
        self.read_type = read_type
        self.out_dir = out_dir
        self.debug = debug
        self.supporting_reads = supporting_reads

    def get_reads(self):

        if self.read_type == 'support':
            bp1_reads, bp1_read_count = self.bp1_supporting_reads()
            bp2_reads, bp2_read_count = self.bp2_supporting_reads()
            bam1 = os.path.join(self.out_dir, "bp1_sv_reads" + ".bam")
            bam2 = os.path.join(self.out_dir, "bp2_sv_reads" + ".bam")
            out = os.path.join(self.out_dir, "sv_support" + ".bam")

            all_reads = bp1_reads + bp2_reads
            total_reads = set(all_reads)

            print("Found %s reads in support of variant" % len(total_reads))

        else:
            bp1_reads, bp1_read_count = self.bp_1_opposing_reads()
            bp2_reads, bp2_read_count = self.bp_2_opposing_reads()
            bam1 = os.path.join(self.out_dir, "bp1_opposing_reads" + ".bam")
            bam2 = os.path.join(self.out_dir, "bp2_opposing_reads" + ".bam")
            out = os.path.join(self.out_dir, "sv_oppose" + ".bam")

            all_reads = bp1_reads + bp2_reads
            total_reads = set(all_reads)

            print("Found %s reads opposing variant" % len(total_reads))

        #-------------------
        # Merge bp1/bp2 reads
        #-------------------
        to_merge = [bam1, bam2]
        merge_bams(out, to_merge)

        return(bp1_reads, bp1_read_count, bp2_reads, bp2_read_count, total_reads)


    def bp1_supporting_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        start = self.bp1 - self.slop
        bp1_reads = []
        out_file = os.path.join(self.out_dir, "bp1_sv_reads" + ".bam")

        with pysam.AlignmentFile(out_file, "wb", template=samfile) as bp1_sv_reads:
            for read in samfile.fetch(self.chrom, start, self.bp1):
                read_end_pos = read.reference_start + read.reference_length
                mate_end_pos = read.next_reference_start + read.reference_length

                if read.is_duplicate:
                    continue

                if not read.is_proper_pair and not read.is_reverse and read.next_reference_start > self.bp2:
                    if self.debug:
                        print("* bp1 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp1_sv_reads.write(read)
                    bp1_reads.append(read.query_name)

                elif read_end_pos == self.bp1:
                    if self.debug:
                        print("* bp1 clipped_read : %s %s [r0: %s, rend: %s]") % (read.query_name, read.seq, read.reference_start, read_end_pos)
                    bp1_sv_reads.write(read)
                    bp1_reads.append(read.query_name)

        count = len(set(bp1_reads))
        pysam.index(out_file)
        return(bp1_reads, count)


    def bp2_supporting_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        end=self.bp2+self.slop
        bp2_reads = []
        out_file = os.path.join(self.out_dir, "bp2_sv_reads" + ".bam")

        with pysam.AlignmentFile(out_file, "wb", template=samfile) as bp2_sv_reads:
            for read in samfile.fetch(self.chrom, self.bp2, end):
                read_end_pos = read.reference_start + read.reference_length
                mate_end_pos = read.next_reference_start + read.reference_length

                if read.is_duplicate:
                    continue

                if not read.is_proper_pair and read.is_reverse and mate_end_pos < self.bp1:
                    if self.debug:
                        print("* bp2 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp2_sv_reads.write(read)
                    bp2_reads.append(read.query_name)

                elif read.reference_start +1 == self.bp2:
                    if self.debug:
                        print("* bp2 clipped_read : %s %s [r0: %s, rend: %s]") % (read.query_name, read.seq, read.reference_start, read_end_pos)
                    bp2_reads.append(read.query_name)
                    bp2_sv_reads.write(read)

        count = len(set(bp2_reads))
        pysam.index(out_file)
        return(bp2_reads, count)


    def bp_1_opposing_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        start = self.bp1 - self.slop
        bp1_reads = []
        # read_names = set()
        out_file = os.path.join(self.out_dir, "bp1_opposing_reads" + ".bam")
        with pysam.AlignmentFile(out_file, "wb", template=samfile) as bp1_opposing_reads:
            for read in samfile.fetch(self.chrom, start, self.bp1):
                read_end_pos = read.reference_start + read.reference_length
                mate_end_pos = read.next_reference_start + read.reference_length

                if read.is_duplicate:
                    continue
                if read.qname in self.supporting_reads:
                    continue

                if (read_end_pos < self.bp1 and read.next_reference_start > self.bp1 and read.next_reference_start < self.bp2) or (read.reference_start > self.bp1 and mate_end_pos < self.bp1):
                    if self.debug:
                        print("* bp1 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp1_opposing_reads.write(read)
                    bp1_reads.append(read.query_name)

                elif read.reference_start < self.bp1 and read_end_pos > self.bp1:
                    if self.debug:
                        print("* bp1 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp1_opposing_reads.write(read)
                    bp1_reads.append(read.query_name)

        count = len(set(bp1_reads))
        pysam.index(out_file)
        return(bp1_reads, count)


    def bp_2_opposing_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        end = self.bp2 + self.slop
        bp2_reads = []
        read_names = set()
        out_file = os.path.join(self.out_dir, "bp2_opposing_reads" + ".bam")
        with pysam.AlignmentFile(out_file, "wb", template=samfile) as bp2_opposing_reads:
            for read in samfile.fetch(self.chrom, self.bp2, end):
                read_end_pos = read.reference_start + read.reference_length
                mate_end_pos = read.next_reference_start + read.reference_length

                if read.is_duplicate:
                    continue
                if read.qname in self.supporting_reads:
                    continue

                if (read_end_pos < self.bp2 and read.next_reference_start > self.bp2) or (read.reference_start > self.bp2 and mate_end_pos < self.bp2 and mate_end_pos > self.bp1):
                    if self.debug:
                        print("* bp2 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp2_opposing_reads.write(read)
                    bp2_reads.append(read.query_name)

                elif read.reference_start < self.bp2 and read_end_pos > self.bp2:
                    if self.debug:
                        print("* bp2 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.query_name, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp2_opposing_reads.write(read)
                    bp2_reads.append(read.query_name)

        count = len(set(bp2_reads))
        pysam.index(out_file)
        return(bp2_reads, count)
