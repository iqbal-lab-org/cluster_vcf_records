import logging
import operator

import pyfastaq

from cluster_vcf_records import vcf_file_read, vcf_record_cluster

class Error (Exception): pass

class VcfClusterer:
    '''Class to cluster one (or more) VCF files. Records in the VCf files
    can be in any order.

    Required parameters:
        vcf_files: list of VCF files
        reference_fasta: FASTA file of reference genome. Must be the same one
           that was used to make the VCF files
        vcf_outfile: name of output VCF file

    Optional parameters:
        max_distance_between_variants: Any variants that are this distance or
           closer will be put in the same cluster. The default of 1 means
           that only adjacent (or overlapping) variants are clustered
        homozygous_only: Set this to True to only load homozygous variants from
           the input VCF files, ie where the genotype is 1/1.
        max_REF_len: When loading the VCF files, any records with REF longer
           than this will be ignored. By default, there is no maximum.
        max_snps_per_cluster: maximum allowed SNPs in one cluster. If a cluster
           has more than this, then it is ignored. Use this option in combination
           with max_REF_len to prevent clusters with huge numbers of SNPs
        source: this is put into the source=foo part of the header of the
           output VCF file.

    Use it like this:

    clusterer = VcfCluster(vcf_files, reference_fasta, vcf_outfile[, options...])
    clusterer.run()
        '''
    def __init__(self, vcf_files, reference_fasta, vcf_outfile,
            max_distance_between_variants=1, homozygous_only=False, max_REF_len=None,
            max_snps_per_cluster=None, source='cluster_vcf_records', merge_method='gramtools'):
        self.vcf_files = vcf_files
        self.reference_seqs = {}
        pyfastaq.tasks.file_to_dict(reference_fasta, self.reference_seqs)
        for seq in self.reference_seqs.values():
            seq.seq = seq.seq.upper()
        self.vcf_outfile = vcf_outfile
        self.max_distance_between_variants = max_distance_between_variants
        self.homozygous_only = homozygous_only
        self.source = source
        self.max_REF_len = max_REF_len
        self.max_snps_per_cluster = max_snps_per_cluster
        self.merge_method = merge_method

        allowed_merge_methods = {'gramtools', 'simple'}
        if self.merge_method not in {'gramtools', 'simple'}:
            raise Error('Erro! merge_method "' + self.merge_method + '" not allowed. Must be one of: ' + ','.join(sorted(list(allowed_merge_methods))))


    @classmethod
    def _load_vcf_files(cls, filename_list, homozygous_only=False, max_REF_len=None):
        '''Loads all the vcf files from filename_list. Returns tuple of:
        1. Sample name. If more than one sample name found, uses the first one
        and warns to stderr
        2. Dictionary. filename => list of header lines for that file
        3. Dictionary. ref name => list of VcfRecords sorted by position'''
        headers = {}
        vcf_records = None
        sample_name = None

        for filename in filename_list:
            headers[filename], new_records = vcf_file_read.vcf_file_to_dict(filename, homozygous_only=homozygous_only, remove_asterisk_alts=True, max_REF_len=max_REF_len, remove_useless_start_nucleotides=True)

            new_sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(headers[filename])
            if sample_name is None and new_sample_name is not None:
                sample_name = new_sample_name
            elif new_sample_name != sample_name:
                logging.warning('Using first sample name found "' + str(sample_name) + '". Found a different (or no) sample name "' + str(new_sample_name) + '", which will not be used')

            if vcf_records is None:
                vcf_records = new_records
            else:
                for ref_name, record_list in new_records.items():
                    if ref_name not in vcf_records:
                        vcf_records[ref_name] = record_list
                    else:
                        vcf_records[ref_name].extend(record_list)

        for record_list in vcf_records.values():
            record_list.sort(key=operator.attrgetter('POS'))

        if sample_name is None:
            logging.warning('No sample name found in VCF files. Going to use "sample"')
            sample_name = 'sample'

        return sample_name, headers, vcf_records


    @classmethod
    def _expand_alts_in_vcf_record_list(cls, vcf_records):
        '''Input: list of vcf_records. Returns new list, where
        any records with >ALT is replaced with one vcf record per ALT.
        This doesn't change FORMAT or INFO columns, which means they
        are now broken for those records'''
        new_vcf_records = []
        for record in vcf_records:
            new_vcf_records.extend(record.to_record_per_alt())
        return new_vcf_records


    @classmethod
    def _cluster_vcf_record_list(cls, vcf_records, max_distance_between_variants=1):
        new_list = [vcf_record_cluster.VcfRecordCluster(max_distance_between_variants=max_distance_between_variants)]

        for vcf_record in vcf_records:
            if not new_list[-1].add_vcf_record(vcf_record):
                new_list.append(vcf_record_cluster.VcfRecordCluster(vcf_record=vcf_record, max_distance_between_variants=max_distance_between_variants))

        return new_list


    def run(self):
        sample_name, vcf_headers, vcf_records = VcfClusterer._load_vcf_files(self.vcf_files, homozygous_only=self.homozygous_only, max_REF_len=self.max_REF_len)

        f_out = pyfastaq.utils.open_file_write(self.vcf_outfile)
        print('##fileformat=VCFv4.2', file=f_out)
        print('##source=', self.source, sep='', file=f_out)
        print('#CHROM', 'POS', 'ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name, sep='\t', file=f_out)

        for ref_name in sorted(vcf_records):
            ref_seq = self.reference_seqs[ref_name]
            cluster_list = VcfClusterer._cluster_vcf_record_list(vcf_records[ref_name], max_distance_between_variants=self.max_distance_between_variants)

            for cluster in cluster_list:
                if self.merge_method == 'gramtools':
                    clustered_vcf = cluster.make_one_merged_vcf_record_for_gramtools(ref_seq, max_snps=self.max_snps_per_cluster)
                    if clustered_vcf is not None:
                        print(clustered_vcf, file=f_out)
                elif self.merge_method == 'simple':
                    clustered_vcf = cluster.make_simple_merged_vcf_with_no_combinations(ref_seq)
                    for vcf in cluster.vcf_records:
                        print(vcf, file=f_out)
                else:
                    raise Error('merge_method "' + self.merge_method + '" not recognised. Cannot continue')

        pyfastaq.utils.close(f_out)

