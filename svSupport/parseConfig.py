import os
import pandas as pd
from worker import worker
from merge_bams import merge_bams
import ntpath


def parse_config(options):
    print("\nExtracting arguments from config file: %s" % options.config)
    base_name = ntpath.basename(options.config)

    # base_name = (os.path.splitext()[0])
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
        options.find_bps = True
        options.guess = df.loc[i, 'guess']

        if df.loc[i, 'chromosome1'] != df.loc[i, 'chromosome2']:
            options.region = df.loc[i, 'chromosome1'] + ":" + str(df.loc[i, 'bp1']) + "-" + df.loc[i, 'chromosome2'] + ":" + str(df.loc[i, 'bp2'])
        else:
            options.region = df.loc[i, 'position']

        bp1, bp2, af, sv_type, integration = worker(options)

        if integration:
            intstring = filter(None, integration)
            intstring = '; '.join(intstring)
            if df.loc[i, 'notes']:
                df.loc[i, 'notes'] = intstring + "; " + df.loc[i, 'notes']
            else: df.loc[i, 'notes'] = intstring

        df.loc[i, 'alf2'] = af
        df.loc[i, 'bp1_c'] = bp1
        df.loc[i, 'bp2_c'] = bp2
        df.loc[i, 'SVtpye'] = sv_type

    mergeAll(options, sample)

    df = df.drop(['bam', 'normal_bam', 'tumour_purity', 'guess', 'sample'], axis=1)
    df.to_csv(outfile, sep="\t", index=False)


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

    allsup = os.path.join(options.out_dir, sample + '_supporting.bam')
    allop = os.path.join(options.out_dir, sample + '_opposing.bam')
    allregions = os.path.join(options.out_dir, sample + '_regions.bam')

    merge_bams(allsup, options.out_dir, su)
    merge_bams(allop, options.out_dir, op)
    merge_bams(allregions, options.out_dir, reg)

