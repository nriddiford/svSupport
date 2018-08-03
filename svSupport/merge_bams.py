import pysam
import os
import ntpath

def merge_bams(out_file, out_dir, bams):
    s_bams = []
    for bam_file in bams:
        sorted_bam = sort_bam(out_dir, bam_file)
        s_bams.append(sorted_bam)

    in_files = ', '.join(s_bams)
    print("Merging bam files %s into '%s'") % (in_files, out_file)
    merge_parameters = ['-f', out_file] + s_bams
    pysam.merge(*merge_parameters)

    sorted_bam = sort_bam(out_dir, out_file)

    # Remove original bams
    for sbam in s_bams:
        try:
            os.remove(sbam)
            os.remove(sbam + ".bai")
        except OSError:
            print("Couldn't remove %s" % sbam)
            pass

    return sorted_bam

def sort_bam(out_dir, bam):
    head, file_name = ntpath.split(bam)
    file_name = os.path.splitext(file_name)[0]
    sorted_bam = os.path.join(out_dir, file_name + ".s" + ".bam")

    try:
        pysam.sort("-o", sorted_bam, bam)
        index_bam(sorted_bam)
    except:
        print("Can't sort %s" % bam)

    try:
        os.remove(bam)
    except:
        print("Can't remove %s" % bam)

    return(sorted_bam)

def index_bam(bam):
    try:
        pysam.index(bam)
    except:
        print("Can't index %s" % bam)
