import os, sys, re
sys.dont_write_bytecode = True
import pandas as pd
from optparse import OptionParser
import ntpath


def makeConfig(options):
    options.outfile = options.sample + '_config.txt'

    with open(options.variants, 'r') as variants:
        purity = get_purity(options)
        sample, group, bamgroup, t_id = getGroup(options.sample)
        sample_bam, normal_bam = getbam(options.bam_dir, bamgroup, group, t_id)
        df = pd.read_csv(variants, delimiter="\t", index_col=False)

        if len(df.index) == 0: sys.exit("No variants in file. Exiting")

        if group in ['D050k', 'D050']:
            sex = 'XX'
        elif group == 'D265' and int(t_id) in [01, 03, 05, 11]:
            print("Setting as female for sample %s" % t_id)
            sex = 'XX'
        elif group == 'D106' and int(t_id) < 23:
            sex = 'XX'
        elif group == 'D197' and str(t_id) in ['09', '11', '13', '15']:
            sex = 'XX'
        else:
            sex = 'XY'

        for i in df.index:
            df.loc[i, 'guess'], df.loc[i, 'normal_bam'] = guess(df.loc[i, 'split_reads'], normal_bam)
            df.loc[i, 'sample'] = sample
            df.loc[i, 'bam'] = sample_bam
            df.loc[i, 'tumour_purity'] = purity
            df.loc[i, 'sex'] = sex


        # TODO add option for Notch filtering
        # df = df[(df['chromosome1'] == '2L') | (df['chromosome2'] == '2L')]
        # df = df[(
        #     (df['chromosome1'] == 'X') & (df['bp1'] > 2500000) & (df['bp1'] < 3500000) |
        #     (df['chromosome2'] == 'X') & (df['bp2'] > 2500000) & (df['bp2'] < 3500000)
        # )]

        df.to_csv(options.outfile, sep="\t", index=False)


def get_purity(options):
    with open(options.purity_file, 'r') as purity_file:
        for l in purity_file:
            parts = l.rstrip().split('\t')
            if options.sample == parts[0]:
                return parts[1]
        print("Can't find corresponding purity for %s in %s" % (options.sample, options.purity_file))
        print("Setting sample purity to 1")
        return 1


def getGroup(sample):
    n = re.search(r'R(.*)', sample)
    if n:
        m = re.search(r'(.*)R', sample)
        t_id = n.group(1)
        group = m.group(1)
    else:
        group = 'HUM'
        n = re.search(r'-(\d+)', sample)
        t_id = n.group(1)

    if group == 'A373':
        bamgroup = 'A370'
    elif group == 'A573':
        bamgroup = 'A572'
    elif group == 'A785-A788':
        bamgroup = 'A785'
    elif group == 'D050' and len(t_id.split('-')) == 1 and int(t_id) >= 10:
        bamgroup = 'D050k'
    else:
        bamgroup = group

    return sample, group, bamgroup, t_id


def guess(sr, nbam):
    if sr == '-':
        return '', nbam
    else:
        return 'T',''


def getbam(bam_dir, bamgroup, group, t_id):
    if group == 'HUM':
        n_id = int(t_id) + 2
        normal_bam = group + "-" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "-" + str(t_id) + '.tagged.filt.SC.RG.bam'

    elif bamgroup == 'D050k':
        n_id = int(t_id) + 1
        normal_bam = group + "R" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

    elif t_id in ['41-1', '41-2']:
        t_no, sid = t_id.split('-')
        n_no = int(t_no) + 1
        n_id = '-'.join(map(str,[n_no, sid]))
        normal_bam = group + "R" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

    elif t_id in ['07-1', '07-2']:
        t_no, sid = t_id.split('-')
        n_no = int(t_no) + 1
        n_id = '-'.join(map(str,[n_no, sid]))
        normal_bam = group + "R" + '0' + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'
    elif group == 'D050':
        n_id = int(t_id) + 1
        normal_bam = group + "R" + '0' + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'
    elif group == 'D197' and str(t_id) in ['01', '03', '05', '07']:
         n_id = int(str(t_id)[-1]) + 1
         normal_bam = group + "R" + '0' + str(n_id) + '.tagged.filt.SC.RG.bam'
         sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

    elif group == 'D265' and str(t_id) in ['01', '03', '05', '07']:
         n_id = int(str(t_id)[-1]) + 1
         normal_bam = group + "R" + '0' + str(n_id) + '.tagged.filt.SC.RG.bam'
         sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

    elif group == 'visor':
        n_id = int(str(t_id)[-1]) + 1
        normal_bam = group + "R" + '0' + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'
    else:
        n_id = int(t_id) + 1
        normal_bam = group + "R" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

    print("Normal: %s Tum: %s", (normal_bam, sample_bam))
    sample_bam = os.path.join(bam_dir, bamgroup, sample_bam)
    normal_bam = os.path.join(bam_dir, bamgroup, normal_bam)

    return sample_bam, normal_bam


def get_args():
    parser = OptionParser()

    parser.add_option("-b",
                      "--bam_dir",
                      dest = "bam_dir",
                      action = "store",
                      help = "Directory containing per-sample bam directories")

    parser.add_option("-v",
                      "--variants",
                      dest = "variants",
                      action = "store",
                      help = "Variants file (as produced by svParser)",
                      metavar = "FILE")

    parser.add_option("-p",
                      "--purity_file",
                      dest = "purity_file",
                      action = "store",
                      help = "File containing 'sample [tab] tumour' purity estimates",
                      metavar="FILE")

    parser.add_option("-s",
                      "--sample",
                      dest="sample",
                      action="store",
                      help="Sample name")

    parser.add_option("-o",
                      "--outfile",
                      dest = "outfile",
                      action = "store",
                      help = "File to annotated variants file to")

    parser.set_defaults(bam_dir='/Volumes/perso/Analysis/Bwa',
                        outfile='/'.join([sys.path[1], '/data/config.txt']),
                        purity_file='/'.join([sys.path[1], '/data/tumour_purity.txt'])
                        )


    options, args = parser.parse_args()

    print(options.purity_file)


    if (options.variants is None):
        parser.print_help()
        print
        sys.exit()

    if options.sample is None:
        options.sample = ntpath.basename(options.variants).split('_')[0]
        print("Extracting sample name from in file: %s" % (options.sample))

    return (options, args)


def main():
    options, args = get_args()

    dir = os.path.dirname(options.outfile)

    if not os.path.exists(dir):
        os.makedirs(dir)

    makeConfig(options)


if __name__ == "__main__":
    sys.exit(main())
