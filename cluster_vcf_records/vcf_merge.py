import datetime

from cluster_vcf_records import vcf_file_read
from cluster_vcf_records import __version__ as cluster_vcf_records_version



def _dict_of_vars_to_vcf_file(variants, outfile):
    '''Input is dict made by vcf_file_read.vcf_file_to_dict_of_vars or
    vcf_file_read.vcf_file_to_dict_of_vars. Output is bare-bones VCF
    file (columns empty wherever possible'''
    header_lines = [
        '##fileformat=VCFv4.2',
        '##source=cluster_vcf_records, version ' + cluster_vcf_records_version,
        '##fileDate=' + str(datetime.date.today()),
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    ]

    with open(outfile, 'w') as f:
        print(*header_lines, sep='\n', file=f)

        for ref_name in sorted(variants):
            for pos in sorted(variants[ref_name]):
                for ref_string in sorted(variants[ref_name][pos]):
                    alts = sorted(list(variants[ref_name][pos][ref_string]))
                    print(ref_name, pos + 1, '.', ref_string, ','.join(alts), '.', 'PASS', 'SVTYPE=MERGED', sep='\t', file=f)


def merge_vcf_files(infiles, outfile, threads=1):
    '''infiles: list of input VCF file to be merge.
    outfile: name of output VCF file.
    threads: number of input files to read in parallel'''
    vars_dict = vcf_file_read.vcf_files_to_dict_of_vars(infiles, threads=threads)
    _dict_of_vars_to_vcf_file(vars_dict, outfile)

