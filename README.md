[![Build Status](https://travis-ci.org/iqbal-lab-org/cluster_vcf_records.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/cluster_vcf_records)

# cluster_vcf_records
Python3 package to cluster VCF records.

Used by gramtools and minos - written specifically for those projects and no others.

## Install
    pip3 install wheel cluster_vcf_records

## Usage
```python
import cluster_vcf_records

input_vcf_file_paths = ['spam.vcf', 'eggs.vcf']
cluster_vcf_records.cluster(input_vcf_file_paths,
                            'reference.fasta',
                            'clustered_output.vcf')
```

For further information and a description of optional arguments, see `help(cluster_vcf_records.vcf_clusterer.VcfClusterer)`

## Develop

Run `python3 -m unittest` from the top-level directory to run all unit tests.

# License
MIT
