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

        if df.loc[i, 'notes'] == '-':
            df.loc[i, 'notes'] = ''

        # if filter_lowFC(i, df):
        #     df.loc[i, 'notes'] = "low read depth ratio"
        #     df.loc[i, 'status'] = 'F'

        genotype = df.loc[i, 'genotype']
        if genotype != 'somatic_tumour': continue

        if df.loc[i, 'chromosome1'] != df.loc[i, 'chromosome2']:
            options.region = df.loc[i, 'chromosome1'] + ":" + str(df.loc[i, 'bp1']) + "-" + df.loc[i, 'chromosome2'] + ":" + str(df.loc[i, 'bp2'])
        else:
            options.region = df.loc[i, 'position']

        if options.guess and df.loc[i, 'status'] != 'F':
            options.find_bps = True

        bp1, bp2, af, sv_type, configuration, notes, split_support, disc_support = worker(options)

        if notes:
            nlist = filter(None, notes)
            nstring = '; '.join(nlist)

            if df.loc[i, 'notes']:
                df.loc[i, 'notes'] = nstring + "; " + df.loc[i, 'notes']
            else: df.loc[i, 'notes'] = nstring

            r = re.compile(".*low read support")
            if filter(r.match, nlist):
                df.loc[i, 'status'] = 'F'
            # Now mark as F if missing read signature
            r = re.compile(".*Missing")
            if filter(r.match, nlist):
                df.loc[i, 'status'] = 'F'
            r = re.compile(".*Contamination")
            if filter(r.match, nlist):
                df.loc[i, 'status'] = 'F'

        if not 'zyg' in sv_type and sv_type != '-':
            df.loc[i, 'type'] = sv_type

        if split_support is not None: df.loc[i, 'split_reads'] = split_support
        if disc_support is not None: df.loc[i, 'disc_reads'] = disc_support

        df.loc[i, 'allele_frequency'] = af
        df.loc[i, 'bp1'] = bp1
        df.loc[i, 'bp2'] = bp2
        df.loc[i, 'configuration'] = configuration
        if df.loc[i, 'chromosome1'] != df.loc[i, 'chromosome2']:
            df.loc[i, 'position'] = df.loc[i, 'chromosome1'] + ":" + str(bp1) + " " + df.loc[i, 'chromosome2'] + ":" + str(bp2)
        else:
            df.loc[i, 'position'] = df.loc[i, 'chromosome1'] + ":" + str(bp1) + "-" + str(bp2)

        if af == 0:
            df.loc[i, 'status'] = 'F'

    mergeAll(options, sample)

    df = df.drop(['bam', 'normal_bam', 'tumour_purity', 'guess', 'sample'], axis=1)
    df = df.sort_values(['chromosome1', 'bp1', 'chromosome2', 'bp2'])
    df.to_csv(outfile, sep="\t", index=False)


def filter_lowFC(i, df):
    c1 = df.loc[i, 'chromosome1']
    rdr = df.loc[i, 'log2(cnv)']

    # if c1 in ['X', 'Y'] and abs(rdr) <

    if abs(rdr) < 0.2:
        return True


def mergeAll(options, sample):
    su = []
    op = []
    reg = []

    for file in os.listdir(options.out_dir):
        if file.endswith("supporting.s.bam"):
            su.append(os.path.join(options.out_dir, file))
        elif file.endswith("opposing.s.bam"):
            op.append(os.path.join(options.out_dir, file))
        elif file.endswith("regions.s.bam"):
            reg.append(os.path.join(options.out_dir, file))

    if(len(su)>1):
        allsup = os.path.join(options.out_dir, sample + '_supporting_dirty.bam')
        allop = os.path.join(options.out_dir, sample + '_opposing.bam')
        allregions = os.path.join(options.out_dir, sample + '_regions.bam')

        sumerged = merge_bams(allsup, options.out_dir, su)
        merge_bams(allop, options.out_dir, op)
        merge_bams(allregions, options.out_dir, reg)

        merged_nodups = os.path.join(sample + '_supporting.bam')
        rmDups(sumerged, merged_nodups, options.out_dir)
