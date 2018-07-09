[![Build Status](https://travis-ci.org/iqbal-lab-org/cluster_vcf_records.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/cluster_vcf_records)

# cluster_vcf_records
Python3 package to cluster VCF records. Used by gramtools and minos - written
specifically for those projects and no others.

## Install
    pip3 install wheel cluster_vcf_records

## Synopsis
To cluster one or more VCF files:

```
import cluster_vcf_records
vcf_files = ['spam.vcf', 'eggs.vcf', 'shrubbery.vcf']
clusterer = cluster_vcf_records.vcf_clusterer.VcfClusterer(vcf_files,
   'reference.fasta',  'out.clustered.vcf')
clusterer.run()
```

For more explanation and a description of optional arguments,
see `help(cluster_vcf_records.vcf_clusterer.VcfClusterer)`

# License
MIT
