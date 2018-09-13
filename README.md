# svSupport

This is a tool to calculate allele frequency of structural variants using a .bam file as input and allowing for the adjustment of allele frequencies based on tumour purity.

svSupport is under constant development. Please feel free to [contact me](mailto:nick.riddiford@curie.fr), or [raise an issue](https://github.com/nriddiford/svSupport/issues) if you encounter any problems.

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

`svSupport -t`

* Bam file: ../data/test.bam
* Chrom: 3L
* bp1: 9892365
* bp2: 9894889
* slop: 500
* search_bps: 0
* debug
* Out dir: test/test_out
