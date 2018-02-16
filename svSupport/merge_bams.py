import pysam
import os

def merge_bams(out_file, bams):
    s_bams = []
    for bam_file in bams:
        try:
            name = os.path.splitext(bam_file)[0]
            print(name)
            sorted_bam = name + ".s" + ".bam"
            pysam.sort("-o", sorted_bam, bam_file)
            pysam.index(sorted_bam)
            s_bams.append(sorted_bam)
        except:
            print("Can't index %s" % bam_file)

    in_files = ', '.join(s_bams)
    print(in_files)
    print("Merging bam files %s into '%s'") % (in_files, out_file)
    merge_parameters = ['-f', out_file] + s_bams
    pysam.merge(*merge_parameters)
    #Remove individual bp files
    for bp_file in s_bams:
        try:
            os.remove(bp_file)
            os.remove(bp_file + ".bai")
        except OSError:
            print("Couldn't remove %s" % bp_file)
            pass

    try:
        pysam.index(out_file)
    except pysam.utils.SamtoolsError:
        print("Can't index %s" % out_file)
        pass
