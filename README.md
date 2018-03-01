# svSupport


```
Usage: svSupport.py [options]

Options:
  -h, --help            show this help message and exit
  -i FILE, --in_file=FILE
                        A sorted .bam file containing the reads supporting the
                        structural variant calls
  -n FILE, --normal_bam=FILE
                        A sorted .bam file for the normal sample used for
                        calculating read depth
  -s SLOP, --slop=SLOP  Distance from breakpoint to look for reads [Default:
                        500]
  -p PURITY, --purity=PURITY
                        Tumour purity e.g. 0.75 [Default: True]
  -f, --find_bps        Look for bps if position not exact [Default: False]
  -l REGION, --loci=REGION
                        The chromosome and breakpoints for a structural
                        variant in the format: 'chrom:bp_1-bp_2'
  -o OUT_DIR, --out_dir=OUT_DIR
                        Directory to write output to [Default: '../out']
  -d, --debug           Run in debug mode [Default: False]
  -t, --test            Run on test data
  -c CONFIG, --config=CONFIG
                        Config file for batch processing sample
                        chromosome:bp1-bp2      purity  ype
  -v VARIANTS_OUT, --variants=VARIANTS_OUT
                        File to write parsed values to
  -g, --guess           Guess type of SV for read searching
```

### Installation
* Create clean environment with [conda](https://conda.io/docs/)
```
pip install conda
conda create -n svSupport python==2.7.12
conda activate svSupport
```
* Install from github
```
git clone https://github.com/nriddiford/svSupport.git
cd svSupport
pip install .
```

### Run on test data

`python svSupport.py -t`

* Bam file: ../data/test.bam
* Chrom: 3L
* bp1: 9892365
* bp2: 9894889
* slop: 500
* search_bps: 0
* debug: True
* test: True
* Out dir: ../test_out
