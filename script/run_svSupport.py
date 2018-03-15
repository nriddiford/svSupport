import fnmatch
import os, sys
sys.dont_write_bytecode = True
import pandas as pd
import re

svs_file='/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/svSupport/data/all_samples_2.txt'

out_file = open('config.txt', 'w')
headers = ['sample', 'event', 'type', 'bam', 'normal_bam', 'locus', 'purity', 'length(Kb)', 'bp1_locus', 'bp2_locus', 'affected_genes', 'guess']

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
	bam_data = os.path.join(bams, group_n)
        depth_data = os.path.join(ratio_files, group_n)
        files = os.listdir(bam_data)
        depth_files = os.listdir(depth_data)
        
#        print("Searcing for ids %s\n%s\n") % (bam_search, normal_sample)
        for name in files:
            if fnmatch.fnmatch(name, bam_search):
                bam_in = os.path.join(bam_data, name)
            elif fnmatch.fnmatch(name, normal_sample):
                normal_bam = os.path.join(bam_data, name)

#	for name in depth_files:
#    	    if fnmatch.fnmatch(name, ratio_search):
#                ratio_file = os.path.join(depth_data, name)

#        print("Found Tumour: %s\nFound Normal: %s\n") % (bam_in, normal_bam)
 
        if "cnv" in variant['source'].lower():
            out_line = [variant['sample'], variant['event'], variant['type'], bam_in, normal_bam, variant['position'], variant['purity'], variant['length(Kb)'], variant['bp1_locus'], variant['bp2_locus'], variant['affected_genes']]
        elif variant['type'] != 'DEL':
            out_line = [variant['sample'], variant['event'], variant['type'], bam_in, 'NA', variant['position'], variant['purity'], variant['length(Kb)'], variant['bp1_locus'], variant['bp2_locus'], variant['affected_genes'], 'T']
        else:
            out_line = [variant['sample'], variant['event'], variant['type'], bam_in, 'NA', variant['position'], variant['purity'], variant['length(Kb)'], variant['bp1_locus'], variant['bp2_locus'], variant['affected_genes']]

        out_file.write('\t'.join(map(str, out_line)) + '\n')
