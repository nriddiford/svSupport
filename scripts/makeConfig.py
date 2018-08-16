import os, sys, re
sys.dont_write_bytecode = True
import pandas as pd
from optparse import OptionParser


def get_purity(options):
    with open(options.purity_file, 'r') as purity_file:
        found = 0
        for l in purity_file:
            parts = l.rstrip().split('\t')
            if options.sample == parts[0]:
                found = 1
                return parts[1]
    if not found:
        print("Can't find corresponding purity for %s in %s" % (options.sample, options.purity_file))
        print("Setting sample purity to 1")
        return 1


def makeConfig(options):
    options.outfile = options.sample + '_config.txt'

    with open(options.variants, 'r') as variants:
        purity = get_purity(options)
        sample, group, bamgroup, t_id = getGroup(options.sample)
        sample_bam, normal_bam = getbam(options.bam_dir, sample, bamgroup, group, t_id)

        df = pd.read_csv(variants, delimiter="\t")

        for i in df.index:
            df.loc[i, 'guess'], df.loc[i, 'normal_bam']  = guess(df.loc[i, 'split_reads'], df.loc[i, 'type'], normal_bam)
            df.loc[i, 'sample'] = sample
            df.loc[i, 'bam'] = sample_bam
            df.loc[i, 'tumour_purity'] = purity

        df.to_csv(options.outfile, sep="\t", index=False)


def guess(sr, t, nbam):
    if sr == '-':
        return '', nbam
    elif t != 'DEL':
        return 'T', ''
    else:
        return '',''


def getGroup(sample):
    n = re.search(r'R(.*)', sample)
    if n:
        m = re.search(r'(.*)R', sample)
        t_id = n.group(1)
        group = m.group(1)
    else:
        group = 'HUM'
        n = re.search(r'-(.*)', sample)
        t_id = n.group(1)

    if group == 'A373':
        bamgroup = 'A370'

    elif group == 'A573':
        bamgroup = 'A572'

    elif group == 'A785-A788':
        bamgroup = 'A785'

    else:
        bamgroup = group

    return sample, group, bamgroup, t_id


def getbam(bam_dir, sample, bamgroup, group, t_id):
    if group == 'HUM':
        n_id = int(t_id) + 2
        normal_bam = group + "-" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "-" + str(t_id) + '.tagged.filt.SC.RG.bam'
    elif sample == 'B241R41-1':
        sample_bam = 'B241R41-1.tagged.filt.SC.RG.bam'
        normal_bam = 'B241R42-1.tagged.filt.SC.RG.bam'
    elif sample == 'B241R41-2':
        sample_bam = 'B241R41-2.tagged.filt.SC.RG.bam'
        normal_bam = 'B241R42-2.tagged.filt.SC.RG.bam'
    else:
        n_id = int(t_id) + 1
        normal_bam = group + "R" + str(n_id) + '.tagged.filt.SC.RG.bam'
        sample_bam = group + "R" + str(t_id) + '.tagged.filt.SC.RG.bam'

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

    parser.set_defaults(bam_dir = '/Users/Nick_curie/Local_data/bam',
                        outfile = 'data/config.txt',
                        variants = '/Users/Nick_curie/Desktop/parserTest/filtered2/summary/merged/HUM-1_annotated_SVs.txt',
                        purity_file = '/Users/Nick_curie/Desktop/script_test/svSupport/data/tumour_purity.txt'
                        )

    options, args = parser.parse_args()

    if (options.variants is None or options.sample is None):
        parser.print_help()
        print
        sys.exit()

    return (options, args)


def main():
    options, args = get_args()

    dir = os.path.dirname(options.outfile)

    if not os.path.exists(dir):
        os.makedirs(dir)

    makeConfig(options)


if __name__ == "__main__":
    sys.exit(main())
