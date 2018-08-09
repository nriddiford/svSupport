import pysam
from collections import defaultdict
import re, os


def disc_reads(read, bp2, forward_reads, reverse_reads):
    if read.is_read1:
        forward_reads[read.query_name] += 1
        if forward_reads[read.query_name] > 1:
            return False

    elif read.is_read2:
        reverse_reads[read.query_name] += 1
        if reverse_reads[read.query_name] > 1:
            return False

    if not read.is_proper_pair and (abs(read.next_reference_start - bp2) <= 350):
        return True


def get_reads(bp_regions, slop, bp_number, chrom, bp, bp2, options, forward_reads, reverse_reads):
    clipped_out = os.path.join(options.out_dir, bp_number + "_clipped_reads" + ".bam")
    disc_out = os.path.join(options.out_dir, bp_number + "_disc_reads" + ".bam")

    samfile = pysam.Samfile(bp_regions, "rb")

    with pysam.AlignmentFile(clipped_out, "wb", template=samfile) as bpReads, pysam.AlignmentFile(disc_out, "wb", template=samfile) as discReads:
        start = bp - slop
        stop = bp + slop

        print ("Looking for reads in %s region: %s:%s-%s" % (bp_number, chrom, start, stop))
        split_reads = 0

        readSig = defaultdict(int)

        for read in samfile.fetch(chrom, start, stop):
            if read.is_reverse:
                direction = 'r'
            else:
                direction = 'f'

            readSig, split_reads, bpID = getClipped(read, bp, direction, bp_number, readSig, split_reads, bpReads)

            if not bpID and disc_reads(read, bp2, forward_reads, reverse_reads):
                tag = ' '.join(['discordant', direction, 'read'])
                tagRead(read, tag)

                if read.is_reverse:
                    bpID = '_'.join([bp_number, direction])
                    print("<---- discordant read supports breakpoints %s: %s") % (bp_number, read.query_name)
                else:
                    print("----> discordant read supports breakpoints %s: %s") % (bp_number, read.query_name)
                    bpID = '_'.join([direction, bp_number])

                readSig[bpID] += 1
                discReads.write(read)


    pysam.index(clipped_out)
    pysam.index(disc_out)
    readSig[None] = 0
    # get most common breakpoint class
    maxValKey = max(readSig, key=readSig.get)

    return (split_reads, maxValKey, clipped_out, disc_out)


def getClipped(read, bp, direction, bp_number, readSig, split_reads, bpReads):
    read_end = read.reference_start + read.reference_length
    bpID = None
    # if clipped
    if re.findall(r'(\d+)[S|H]', read.cigarstring):
        # if right-clipped
        if bp == read_end:
            bpID = rightClipped(read, direction, bp_number)
        # if read is left-clipped
        elif not bpID and bp == read.reference_start + 1:
            bpID = leftClipped(read, direction, bp_number)

    if bpID:
        tag = ' '.join(['clipped', direction, 'read'])
        tagRead(read, tag)
        readSig[bpID] += 1
        split_reads += 1
        bpReads.write(read)

    return readSig, split_reads, bpID


def tagRead(read, tag):
    read.set_tag('SV', tag, value_type='Z')
    return read


def rightClipped(read, direction, bp_number):
    if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
        if direction == 'f':
            print("---> read clipped to right: %s --[-> %s") % (read.cigarstring, read.query_name)
        else:
            print("<--- read clipped to right: %s <--[- %s") % (read.cigarstring, read.query_name)
        bpID = '_'.join([direction, bp_number])
        return bpID


def leftClipped(read, direction, bp_number):
    if re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
        if direction == 'f':
            print("---> read clipped to left: %s -]--> %s") % (read.cigarstring, read.query_name)
        else:
            print("<--- read clipped to left: %s <-]-- %s") % (read.cigarstring, read.query_name)
        bpID = '_'.join([bp_number, direction])
        return bpID