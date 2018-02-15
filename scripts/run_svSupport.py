import fnmatch
import os, sys
sys.dont_write_bytecode = True
import pandas as pd

import re

svs_file='/Users/Nick_curie/Desktop/svParser/filtered/summary/merged/all_samples.txt'

out_file = open('config.txt', 'w')
headers = ['sample', 'locus', 'purity', 'read_depth']

out_file.write('\t'.join(headers) + '\n')
data = '/Users/Nick_curie/Desktop/script_test/svSupport/data'
tumour_purity = '/Users/Nick_curie/Desktop/tumour_purity.txt'


def get_purity():
    with open(tumour_purity, 'r') as purity_file:
        sample_purity = {}
        for l in purity_file:
            parts = l.rstrip().split('\t')
            sample_purity[parts[0]] = parts[1]
    return(sample_purity)


print(headers)
with open(svs_file, 'r') as variants:
    dataset=pd.read_csv(variants,delimiter="\t")

    purity = get_purity()

    for index, variant in dataset.iterrows():
        if variant['genotype'] != 'somatic_tumour':
            continue
        if variant['type'] != 'DEL':
            continue

        ratio_search = variant['sample'] + '*_ratio.txt'
        bam_search = variant['sample'] + '*.bam'

        sample = variant['sample']
        m = re.search(r'(.*)R',sample)

        group = m.group(1)

        print(group)


        if sample in purity:
            variant['purity'] = purity[sample]
        else:
            print("Can't find corresponding purity for %s" % sample)
            variant['purity'] = 1

        ratio_file = None
        bam_in = None
        files = os.listdir(data)

        for name in files:
            if fnmatch.fnmatch(name, ratio_search):
                ratio_file = os.path.join(data, name)
            if fnmatch.fnmatch(name, bam_search):
                bam_in = os.path.join(data, name)

        if "cnv" in variant['source'].lower():
            print("CNV-Seq var: %s, %s, %s, %s") % (variant['sample'], ratio_file, bam_in, variant['purity'])
            print(variant['sample'], bam_in, variant['position'], variant['purity'], ratio_file)
        else:
            print(variant['sample'], bam_in, variant['position'], variant['purity'], 'NA')
