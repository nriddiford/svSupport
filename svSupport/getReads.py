import pysam
from collections import defaultdict
import re, os


def get_reads(bp_regions, bp_number, chrom, bp, bp2, options, forward_reads, reverse_reads, chroms, bp_sig):
    """Get reads supporting a breakpoint and write both discordant and clipped supporting reads.
       Count the evidence supporting different types of breakpoint arrangements and return the one
       with most support:

       f_bp : ---->[x]=====
       bp_f : =====[x]---->
       r_bp : <----[x]=====
       bp_r : =====[x]<----

       """
    clipped_out = os.path.join(options.out_dir, bp_number + "_clipped_reads" + ".bam")
    disc_out = os.path.join(options.out_dir, bp_number + "_disc_reads" + ".bam")

    samfile = pysam.Samfile(bp_regions, "rb")

    with pysam.AlignmentFile(clipped_out, "wb", template=samfile) as bpReads, pysam.AlignmentFile(disc_out, "wb", template=samfile) as discReads:
        split_reads = 0

        start = bp - 500
        stop = bp + 500

        alien_split = 0
        alien_paired = 0
        te_tagged = defaultdict(int)

        for read in samfile.fetch(chrom, start, stop):
            if not read.infer_read_length():
                """This is a problem - and will skip over reads with no mapped mate (which also have no cigar)"""
                continue
            if read.is_duplicate:
                continue

            if read.is_reverse:
                direction = 'r'
            else:
                direction = 'f'

            read, bp_sig, split_reads, bpID = getClipped(read, bp, direction, bp_number, bp_sig, split_reads, options)

            if bpID:
                bpReads.write(read)

            if disc_reads(read, bp2, forward_reads, reverse_reads):
                tag = ' '.join(['discordant', direction, 'read'])
                tagRead(read, tag)
                if read.is_reverse:
                    bpID = '_'.join([bp_number, direction])
                    if options.debug:
                        print("<---- discordant read supports breakpoints %s: %s") % (bp_number, read.query_name)
                else:
                    bpID = '_'.join([direction, bp_number])
                    if options.debug:
                        print("----> discordant read supports breakpoints %s: %s") % (bp_number, read.query_name)


                bp_sig[bpID] += 1
                discReads.write(read)
            if options.chromfile:
                alien_split, alien_paired = getAlienDNA(bpID, read, alien_split, alien_paired, chroms, bp_number, direction)
            te_tagged = getTaggedReads(bpID, read, te_tagged, bp_number, direction)



    pysam.index(clipped_out)
    pysam.index(disc_out)
    bp_sig[None] = 0
    te_tagged[None] = 0

    # get most common breakpoint class
    maxValKey = max(bp_sig, key=bp_sig.get)
    te = max(te_tagged, key=te_tagged.get)

    return (split_reads, maxValKey, clipped_out, disc_out, alien_split, alien_paired, te, te_tagged, bp_sig)


def getAlienDNA(bpID, read, alien_split, alien_paired, chroms, bp_number, direction):
    if bpID and read.is_supplementary:
        try:
            sa_chrom = read.get_tag('SA').split(',')[0]
            if sa_chrom not in chroms:
                # print "Clipped reads partially aligns to non-reference chromosome: %s (%s)" % (sa_chrom, read.query_name)
                alien_split += 1
        except:
            pass

    elif read.next_reference_name not in chroms:
        if bp_number == 'bp1' and direction == 'f':
            print "forward read has mate mapped to non-reference chromosome: %s (%s)" % (read.next_reference_name, read.query_name)
            alien_paired += 1
        elif bp_number == 'bp2' and direction == 'r':
            print "reverse read has mate mapped to non-reference chromosome: %s (%s)" % (read.next_reference_name, read.query_name)
            alien_paired += 1

    return alien_split, alien_paired


def getTaggedReads(bpID, read, te_tagged, bp_number, direction):
    if bpID and read.is_supplementary:
        try:
            sa_te = read.get_tag('AD').split(',')[0]
            sa_te = '_'.join(sa_te.split('_')[1:])
            print "clipped:", sa_te, read.query_name
            te_tagged[sa_te] += 1
        except KeyError:
            pass
    try:
        te = read.get_tag('BR').split(',')[0]
        te = '_'.join(te.split('_')[1:])
        if bp_number == 'bp1' and direction == 'f':
            print "foward read has mate in te:", te, read.query_name
            te_tagged[te] += 1
        elif bp_number == 'bp2' and direction == 'r':
            print "reverse read has mate in te:", te, read.query_name
            te_tagged[te] += 1
    except KeyError:
        pass

    return te_tagged


def getClipped(read, bp, direction, bp_number, bp_sig, split_reads, options):
    """If reads are clipped at the breakpoint, classify them based on their mapping signature:
       f_bp : ---->[x]=====
       bp_f : =====[x]---->
       r_bp : <----[x]=====
       bp_r : =====[x]<----
       """

    bpID = None
    # if clipped
    if re.findall(r'(\d+)[S|H]', read.cigarstring):
        # if right-clipped
        if bp == read.reference_end:
            bpID = rightClipped(read, direction, bp_number, options)
        # if read is left-clipped
        elif not bpID and bp == read.reference_start + 1:
            bpID = leftClipped(read, direction, bp_number, options)

    if bpID:
        tag = ' '.join(['clipped', direction, 'read'])
        tagRead(read, tag)
        bp_sig[bpID] += 1
        split_reads += 1

    return read, bp_sig, split_reads, bpID


def rightClipped(read, direction, bp_number, options):
    """Looks for reads that are clipped to the right of breakpoint"""
    if re.findall(r".*?M(\d+)[S|H]", read.cigarstring):
        if options.debug:
            if direction == 'f':
                print("---> read clipped to right: %s --[-> %s") % (read.cigarstring, read.query_name)
            else:
                print("<--- read clipped to right: %s <--[- %s") % (read.cigarstring, read.query_name)
        bpID = '_'.join([direction, bp_number])
        return bpID


def leftClipped(read, direction, bp_number, options):
    """Looks for reads that are clipped to the left of breakpoint"""
    if re.findall(r'(\d+)[S|H].*?M', read.cigarstring):
        if options.debug:
            if direction == 'f':
                print("---> read clipped to left: %s -]--> %s") % (read.cigarstring, read.query_name)
            else:
                print("<--- read clipped to left: %s <-]-- %s") % (read.cigarstring, read.query_name)
        bpID = '_'.join([bp_number, direction])
        return bpID


def disc_reads(read, bp2, forward_reads, reverse_reads):
    """Return true if read is both discordant (not read.is_proper_pair)
       and not contained in f/r reads dict. This prevents the same read being counted
       twice in regions that overlap"""
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


def tagRead(read, tag):
    """Tag read with custom tag 'SV' - indicating its involvement in SV"""
    read.set_tag('SV', tag, value_type='Z')
    return read
