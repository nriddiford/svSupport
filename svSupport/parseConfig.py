import os, re
import pandas as pd
from worker import worker
from merge_bams import merge_bams
from worker import rmDups
import ntpath


def parse_config(options):
    print("\nExtracting arguments from config file: %s" % options.config)
    base_name = ntpath.basename(options.config)

    if not options.variants_out:
        sample = base_name.split('_')[0]
        outfile = sample + '_svSupport.txt'
        options.variants_out = outfile

    df = pd.read_csv(options.config, delimiter="\t")
    df = df.where((pd.notnull(df)), None)

    for i in df.index:
        options.in_file = df.loc[i, 'bam']
        options.purity = float(df.loc[i, 'tumour_purity'])
        options.normal_bam = df.loc[i, 'normal_bam']
        options.guess = df.loc[i, 'guess']
        options.sex = df.loc[i, 'sex']

        if df.loc[i, 'notes'] == '-': df.loc[i, 'notes'] = ''
        if df.loc[i, 'status'] == '-': df.loc[i, 'status'] = ''

        genotype = df.loc[i, 'genotype']
        if genotype != 'somatic_tumour': continue

        if df.loc[i, 'chromosome1'] != df.loc[i, 'chromosome2']:
            options.region = df.loc[i, 'chromosome1'] + ":" + str(df.loc[i, 'bp1']) + "-" + df.loc[i, 'chromosome2'] + ":" + str(df.loc[i, 'bp2'])
        else:
            options.region = df.loc[i, 'position']

        # TODO this can be cleaned up now (seeing as we're not marking vars prior to svSupport
        if options.guess and df.loc[i, 'status'] != 'F':
            options.find_bps = True

        bp1, bp2, af, sv_type, configuration, notes, split_support, disc_support = worker(options)

        if options.normal_bam:
            df.loc[i, 'configuration'] = sv_type
            sv_type = df.loc[i, 'type']
        else:
             df.loc[i, 'configuration'] = configuration

        notes = mark_low_FC(notes, options.sex, df.loc[i, 'log2(cnv)'], sv_type, df.loc[i, 'chromosome1'], split_support)

        nlist = filter(None, notes)
        nstring = '; '.join(nlist)
        if df.loc[i, 'notes']:
            df.loc[i, 'notes'] = nstring + "; " + df.loc[i, 'notes']
        else: df.loc[i, 'notes'] = nstring

        df.loc[i, 'status'] = mark_filters(notes)
        df.loc[i, 'type'] = sv_type

        if split_support is not None: df.loc[i, 'split_reads'] = split_support
        if disc_support is not None: df.loc[i, 'disc_reads'] = disc_support

        df.loc[i, 'allele_frequency'] = af
        df.loc[i, 'bp1'] = bp1
        df.loc[i, 'bp2'] = bp2

        if df.loc[i, 'chromosome1'] != df.loc[i, 'chromosome2']:
            df.loc[i, 'position'] = df.loc[i, 'chromosome1'] + ":" + str(bp1) + " " + df.loc[i, 'chromosome2'] + ":" + str(bp2)
        else:
            df.loc[i, 'position'] = df.loc[i, 'chromosome1'] + ":" + str(bp1) + "-" + str(bp2)

        if af == 0:
            df.loc[i, 'status'] = 'F'

    mergeAll(options, sample)

    df = df.drop(['bam', 'normal_bam', 'tumour_purity', 'guess', 'sample', 'sex'], axis=1)
    df = df.sort_values(['chromosome1', 'bp1', 'chromosome2', 'bp2'])
    df.to_csv(outfile, sep="\t", index=False)


def mark_low_FC(notes, sex, fc, sv_type, chrom, split_support):
    if split_support >= 5:
        return notes
    elif split_support:
        hom_pass = 0.2
        het_pass = 0.1
    else:
        hom_pass = 0.6
        het_pass = 0.4

    if sv_type in ['DEL', 'DUP', 'TANDUP']:
        if chrom in ['X', 'Y'] and sex == 'XY':
            if abs(fc) < hom_pass:
                notes.append("low FC")
        elif abs(fc) < het_pass:
            notes.append("low FC")

    return notes


def mark_filters(notes):
    filters = ['low read support', 'missing', 'contamination', 'low depth', 'low FC', 'excluded']
    for n in notes:
        for f in filters:
            if f in n:
                return 'F'

def mergeAll(options, sample):
    su = []
    op = []
    reg = []
    print("MergeAll")
    for file in os.listdir(options.out_dir):
        if file.endswith("supporting.s.bam"):
            su.append(os.path.join(options.out_dir, file))
        elif file.endswith("opposing.s.bam"):
            op.append(os.path.join(options.out_dir, file))
        elif file.endswith("regions.s.bam"):
            reg.append(os.path.join(options.out_dir, file))

    if len(su) > 1:
        allsup = os.path.join(options.out_dir, sample + '_supporting_dirty.bam')
        allop = os.path.join(options.out_dir, sample + '_opposing_dirty.bam')
        allregions = os.path.join(options.out_dir, sample + '_regions_dirty.bam')

        sumerged = merge_bams(allsup, options.out_dir, su)
        opmerged = merge_bams(allop, options.out_dir, op)
        remerged = merge_bams(allregions, options.out_dir, reg)

        merged_su_nodups = os.path.join(sample + '_supporting.bam')
        merged_op_nodups = os.path.join(sample + '_opposing.bam')
        merged_regs_nodups = os.path.join(sample + '_regions.bam')

        rmDups(sumerged, merged_su_nodups, options.out_dir)
        rmDups(opmerged, merged_op_nodups, options.out_dir)
        rmDups(remerged, merged_regs_nodups, options.out_dir)

