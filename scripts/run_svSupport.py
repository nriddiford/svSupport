import fnmatch
import os, sys
sys.dont_write_bytecode = True
import pandas as pd
import re
from optparse import OptionParser

def get_purity(options):
    with open(options.purity_file, 'r') as purity_file:
        sample_purity = {}
        for l in purity_file:
            parts = l.rstrip().split('\t')
            if options.sample == parts[0]:
                sample_purity[parts[0]] = parts[1]
                return(sample_purity)


def makeConfig(options):
    with open(options.variants, 'r') as variants, open(options.outfile, 'w') as config_out:

        headers = ['event', 'type', 'bam', 'normal_bam', 'locus', 'purity', 'length(Kb)', 'guess']
        config_out.write('\t'.join(headers) + '\n')

        dataset = pd.read_csv(variants, delimiter="\t")
        purity = get_purity(options)

        for index, variant in dataset.iterrows():
            if variant['genotype'] != 'somatic_tumour':
                continue
            if variant['chromosome1'] != variant['chromosome2']:
                continue

            ratio_search = variant['sample'] + '*_ratio.txt'
            bam_search = variant['sample'] + '.tagged.filt.SC.RG.bam'

            sample = variant['sample']
            n = re.search(r'R(.*)',sample)
            if n:
                R_id = n.group(1)

            m = re.search(r'(.*)R',sample)

            if m:
                group = m.group(1)
            else:
                group = 'HUM'
                n = re.search(r'-(.*)',sample)
                R_id = n.group(1)

            if group == 'A373':
                group_n='A370'

            elif group == 'A573':
                    group_n='A572'

            elif group == 'A785-A788':
                 group_n = 'A785'

            else:
                 group_n = group

            if group == 'HUM':
                R_id = int(R_id) + 2
                normal_sample = group + "-" + str(R_id) + '.tagged.filt.SC.RG.bam'
            elif sample == 'B241R41-1':
                normal_sample = 'B241R42-1.tagged.filt.SC.RG.bam'
            elif sample == 'B241R41-2':
                normal_sample = 'B241R42-2.tagged.filt.SC.RG.bam'
            else:
                R_id = int(R_id) + 1
                normal_sample = group + "R" + str(R_id) + '.tagged.filt.SC.RG.bam'

            if sample in purity:
                variant['purity'] = purity[sample]
            else:
                print("Can't find corresponding purity for %s" % sample)
                print(sample)
                variant['purity'] = 1

            normal_bam = None
            bam_in = None
            bam_data = os.path.join(options.bam_dir, group_n)
            files = os.listdir(bam_data)

            for name in files:
                if fnmatch.fnmatch(name, bam_search):
                    bam_in = os.path.join(bam_data, name)
                elif fnmatch.fnmatch(name, normal_sample):
                    normal_bam = os.path.join(bam_data, name)

            if "cnv" in variant['source'].lower():
                out_line = [variant['event'], variant['type'], bam_in, normal_bam, variant['position'], variant['purity']]
            elif variant['type'] != 'DEL':
                out_line = [variant['event'], variant['type'], bam_in, 'NA', variant['position'], variant['purity'], 'T']
            else:
                out_line = [variant['event'], variant['type'], bam_in, 'NA', variant['position'], variant['purity'] ]

            config_out.write('\t'.join(map(str, out_line)) + '\n')

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

    parser.add_option("-r",
                      "--ratio_dir",
                      dest="ratio_dir",
                      action="store",
                      help="Directory containing per-sample ratio directories")

    parser.add_option("-o",
                      "--outfile",
                      dest = "outfile",
                      action = "store",
                      help = "File to annotated variants file to")

    parser.set_defaults(bam_dir = '/Users/Nick_curie/Local_data/bam',
                        outfile = 'data/config.txt',
                        purity_file = '/Users/Nick_curie/Desktop/script_test/svSupport/data/tumour_purity.txt'
                        )

    options, args = parser.parse_args()

    return (options, args)


def main():
    options, args = get_args()

    dir = os.path.dirname(options.outfile)

    if not os.path.exists(dir):
        os.makedirs(dir)

    if options.variants is None:
        parser.print_help()
        print()
    else:
        makeConfig(options)


if __name__ == "__main__":
    sys.exit(main())
