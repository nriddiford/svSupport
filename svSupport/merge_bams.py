import pysam
import os
import ntpath

def merge_bams(out_file, out_dir, bams):
    s_bams = []
    for bam_file in bams:
        try:
            head, file_name = ntpath.split(bam_file)
            sorted_bam = os.path.join(out_dir, file_name + ".s" + ".bam")
            pysam.sort("-o", sorted_bam, bam_file)
            os.remove(bam_file)
            s_bams.append(sorted_bam)
            pysam.index(sorted_bam)
        except:
            print("Can't sort %s" % bam_file)

    in_files = ', '.join(s_bams)
    print("Merging bam files %s into '%s'") % (in_files, out_file)
    merge_parameters = ['-f', out_file] + s_bams
    pysam.merge(*merge_parameters)

    head, file_name = ntpath.split(out_file)
    merged_bam = os.path.splitext(file_name)[0]

    sorted_bam = os.path.join(out_dir, merged_bam + ".s" + ".bam")
    pysam.sort("-o", sorted_bam, out_file)
    pysam.index(out_file)

    #Remove individual bp files
    for bp_file in s_bams:
        try:
            abs_file = os.path.join(out_dir, bp_file)
            os.remove(abs_file)
            os.remove(abs_file + ".bai")
        except OSError:
            print("Couldn't remove %s" % abs_file)
            pass

    return(sorted_bam)
