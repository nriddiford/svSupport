import pandas as pd
from svSupport import worker


def parse_config(options):
    print("\nExtracting arguments from config file: %s" % options.config)
    base_name = (os.path.splitext(options.config)[0])
    if not options.variants_out:
        sample = base_name.split('_')[0]
        outfile = sample + '_svSupport.txt'
        options.variants_out = outfile

    df = pd.read_csv(options.config, delimiter="\t")
    df = df.where((pd.notnull(df)), None)

    for i in df.index:
        options.in_file = df.loc[i, 'bam']

        options.region = df.loc[i, 'position']
        options.purity = float(df.loc[i, 'tumour_purity'])
        options.normal_bam = df.loc[i, 'normal_bam']
        options.find_bps = True
        options.guess = df.loc[i, 'guess']

        # if df.loc[i, 'type'] in ['DEL', 'DUP']:
        if df.loc[i, 'chromosome1'] == df.loc[i, 'chromosome2']:
            chrom, bp1, bp2, allele_frequency = worker(options)
        else:
            chrom, bp1, bp2, allele_frequency = (0,0,0,0)

        df.loc[i, 'alf2'] = allele_frequency
        df.loc[i, 'bp1_c'] = bp1
        df.loc[i, 'bp2_c'] = bp2

    df = df.drop(['bam', 'normal_bam', 'tumour_purity', 'guess', 'sample'], axis=1)
    df.to_csv(outfile, sep="\t", index=False)