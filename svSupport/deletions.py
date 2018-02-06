import pysam
import os, re
from merge_bams import merge_bams

class Deletions(object):
    def __init__(self, bam_in, chrom, bp1, bp2, slop, read_type, out_dir, debug=False):

        self.bam_in = bam_in
        self.chrom = chrom
        self.bp1 = bp1
        self.bp2 = bp2
        self.slop = slop
        self.read_type = read_type
        self.reads = [self.bam_in, self.chrom, self.bp1, self.bp2, self.slop]
        self.out_dir = out_dir
        self.debug = debug


    def get_reads(self):

        # self.read_type = read_type

        if self.read_type == 'support':
            bp1_reads, bp1_read_count = self.bp1_supporting_reads()
            bp2_reads, bp2_read_count = self.bp2_supporting_reads()
            bam1 = os.path.join(self.out_dir, "bp1_sv_reads" + ".bam")
            bam2 = os.path.join(self.out_dir, "bp2_sv_reads" + ".bam")
            out = os.path.join(self.out_dir, "sv_support" + ".bam")
            total_support = bp1_read_count + bp2_read_count

            print("Found %s reads in support of variant" % total_support)

        else:
            bp1_reads, bp1_read_count = self.bp_1_opposing_reads()
            bp2_reads, bp2_read_count = self.bp_2_opposing_reads()
            bam1 = os.path.join(self.out_dir, "bp1_opposing_reads" + ".bam")
            bam2 = os.path.join(self.out_dir, "bp2_opposing_reads" + ".bam")
            out = os.path.join(self.out_dir, "sv_oppose" + ".bam")
            total_oppose = bp1_read_count + bp2_read_count

            print("Found %s reads opposing variant" % total_oppose)

        to_merge = [bam1, bam2]
        merge_bams(out, to_merge)
        return(bp1_reads, bp1_read_count, bp2_reads, bp2_read_count)


    def bp1_supporting_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        start = self.bp1 - self.slop

        bp1_reads = []
        out_file = os.path.join(self.out_dir, "bp1_sv_reads" + ".bam")

        count = 0
        with pysam.AlignmentFile(out_file, "wb", template=samfile) as bp1_sv_reads:
            for read in samfile.fetch(self.chrom, start, self.bp1):
                read_end_pos = read.reference_start + read.alen
                mate_end_pos = read.next_reference_start + read.alen

                if read.is_duplicate:
                    print(read.qname)
                    continue

                if not read.is_proper_pair and not read.is_reverse and read.next_reference_start > self.bp2:
                    if self.debug:
                        print("* bp1 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                    bp1_sv_reads.write(read)
                    bp1_reads.append(read.qname)
                    count += 1

                elif read_end_pos == self.bp1:
                    if self.debug:
                        print("* bp1 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.reference_start, read_end_pos)
                    bp1_sv_reads.write(read)
                    bp1_reads.append(read.qname)
                    count += 1

        pysam.index(out_file)

        return(bp1_reads, count)


    def bp2_supporting_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        end=self.bp2+self.slop

        bp2_reads = []
        out_file = os.path.join(self.out_dir, "bp2_sv_reads" + ".bam")
        bp2_sv_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)

        count = 0

        for read in samfile.fetch(self.chrom, self.bp2, end):
            read_end_pos = read.reference_start + read.alen
            mate_end_pos = read.next_reference_start + read.alen

            if read.is_duplicate:
                continue

            if not read.is_proper_pair and read.is_reverse and mate_end_pos < self.bp1:
                if self.debug:
                    print("* bp2 disc_read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                bp2_sv_reads.write(read)
                bp2_reads.append(read.qname)
                count += 1

            elif read.reference_start +1 == self.bp2:
                if self.debug:
                    print("* bp2 clipped_read : %s %s [r0: %s, rend: %s]") % (read.qname, read.seq, read.reference_start, read_end_pos)
                bp2_reads.append(read.qname)
                bp2_sv_reads.write(read)
                count += 1

        bp2_sv_reads.close()
        pysam.index(out_file)

        return(bp2_reads, count)


    def bp_1_opposing_reads(self):

        samfile = pysam.Samfile(self.bam_in, "rb")
        start = self.bp1 - self.slop
        bp1_reads = []
        read_names = set()

        out_file = os.path.join(self.out_dir, "bp1_opposing_reads" + ".bam")
        bp1_opposing_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
        count = 0

        for read in samfile.fetch(self.chrom, start, self.bp1):
            read_end_pos = read.reference_start + read.alen
            mate_end_pos = read.next_reference_start + read.alen

            if read.is_duplicate:
                continue

            if read.is_proper_pair and not read.is_reverse and read.next_reference_start > self.bp1 and not read.is_supplementary and not read_end_pos == self.bp1:
                if self.debug:
                    print("* bp1 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                bp1_opposing_reads.write(read)
                bp1_reads.append(read.qname)
                read_names.update(read.qname)
                count += 1

            elif read.reference_start < self.bp1 and read_end_pos > self.bp1:
                if self.debug:
                    print("* bp1 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                bp1_opposing_reads.write(read)
                bp1_reads.append(read.qname)
                read_names.update(read.qname)
                count += 1

      # print(read_names)
        bp1_opposing_reads.close()
        pysam.index(out_file)

        return(bp1_reads, count)


    def bp_2_opposing_reads(self):
        samfile = pysam.Samfile(self.bam_in, "rb")
        end = self.bp2 + self.slop
        # self.bp1 = bp1
        # self.bp2 = bp2
        # samfile = pysam.Samfile(bamFile, "rb")
        # end=bp2+slop
        bp2_reads = []
        out_file = os.path.join(self.out_dir, "bp2_opposing_reads" + ".bam")
        bp2_opposing_reads = pysam.AlignmentFile(out_file, "wb", template=samfile)
        count = 0

        for read in samfile.fetch(self.chrom, self.bp2, end):
            read_end_pos = read.reference_start + read.alen
            mate_end_pos = read.next_reference_start + read.alen

            if read.is_duplicate:
                continue

            if read.is_proper_pair and read.is_reverse and read.next_reference_start < self.bp2 and not read.is_supplementary and read.reference_start +1 != self.bp2 and read.next_reference_start +1 != self.bp2:
                if self.debug:
                    print("* bp2 opposing read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                bp2_opposing_reads.write(read)
                bp2_reads.append(read.qname)
                count += 1

            elif read.reference_start > self.bp2 and read_end_pos < self.bp2:
                if self.debug:
                    print("* bp2 spanning read    : %s %s [rs:e: %s-%s, ms:e: %s-%s]") % (read.qname, read.seq, read.reference_start, read_end_pos, read.next_reference_start, mate_end_pos)
                bp2_opposing_reads.write(read)
                bp2_reads.append(read.qname)
                count += 1

        bp2_opposing_reads.close()
        pysam.index(out_file)

        return(bp2_reads, count)