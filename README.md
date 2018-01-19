# svSupport


```
Usage: svSupport.py [options]

Options:
  -h, --help            show this help message and exit
  -i FILE, --inFile=FILE
                        A sorted .bam file containing the reads supporting the
                        structural variant calls
  -s SLOP, --slop=SLOP  Distance from breakpoint to look for readsDefault: 500
  -l REGION, --loci=REGION
                        The chromosome and breakpoints for a structural
                        variant in the format: 'chrom:bp_1-bp_2'
```

# Run on test data
* Genome: D_mel_6.12
* 2.5 kb DEL
* chrom: 3L
* bp1: 9892365
* bp2:9894889
* allele frequency: 0.36

`python svSupport.py -i ../data/test.bam -l 3L:9892365-9894889 -d -o ../test_out -s 500`
