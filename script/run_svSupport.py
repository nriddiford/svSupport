import fnmatch
import os, sys
sys.dont_write_bytecode = True
import pandas as pd
import re

svs_file='/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/svSupport/data/all_samples_2.txt'

out_file = open('config.txt', 'w')
headers = ['sample', 'event', 'type', 'bam', 'locus', 'purity', 'read_depth', 'length(Kb)', 'bp1_locus', 'bp2_locus', 'affected_genes']

out_file.write('\t'.join(headers) + '\n')
bams = '/data/kdi_prod/project_result/948/01.00/Analysis/Bwa'
tumour_purity = '/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Control-Freec/tumour_purity.txt'
ratio_files = '/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Control-Freec'

def get_purity():
    with open(tumour_purity, 'r') as purity_file:
        sample_purity = {}
        for l in purity_file:
            parts = l.rstrip().split('\t')
            sample_purity[parts[0]] = parts[1]
    return(sample_purity)


with open(svs_file, 'r') as variants:
    dataset=pd.read_csv(variants,delimiter="\t")

    purity = get_purity()

    for index, variant in dataset.iterrows():
        if variant['genotype'] != 'somatic_tumour':
            continue
        if variant['type'] != 'DEL':
            continue
	
        ratio_search = variant['sample'] + '*_ratio.txt'
        bam_search = variant['sample'] + '.tagged.filt.SC.RG.bam'

        sample = variant['sample']
	   
        m = re.search(r'(.*)R',sample)
	
        if m:
            group = m.group(1)
	else:
            group = 'HUM'

	if group == 'A373':
	    group='A370'
	
	elif group == 'A573':
            group='A572'
  
	elif group == 'A785-A788':
	     group = 'A785'

        if sample in purity:
            variant['purity'] = purity[sample]
        else:
            print("Can't find corresponding purity for %s" % sample)
            print(sample)
            variant['purity'] = 1

        ratio_file = None
        bam_in = None
	bam_data = os.path.join(bams, group)
        depth_data = os.path.join(ratio_files, group)
        files = os.listdir(bam_data)
        depth_files = os.listdir(depth_data)

        for name in files: 
            if fnmatch.fnmatch(name, bam_search):
                bam_in = os.path.join(bam_data, name)

	for name in depth_files:
    	    if fnmatch.fnmatch(name, ratio_search):
                ratio_file = os.path.join(depth_data, name)
        
        if "cnv" in variant['source'].lower():
            out_line = [variant['sample'], variant['event'], variant['type'],  bam_in, variant['position'], variant['purity'], ratio_file, variant['length(Kb)'], variant['bp1_locus'], variant['bp2_locus'], variant['affected_genes']]
        else:
            out_line = [variant['sample'], variant['event'], variant['type'], bam_in, variant['position'], variant['purity'], 'NA', variant['length(Kb)'], variant['bp1_locus'], variant['bp2_locus'], variant['affected_genes']]

        out_file.write('\t'.join(map(str, out_line)) + '\n')
        
