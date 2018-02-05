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
        min_SNP_qual: Set this to the minimum QUAL score for inclusion of variant call from samtools
        min_dp4: Set this to the minimum high quality read depth on both strands for inclusion of a variant call from samtools
        min_GT_conf: Set this to the minimum genotype confidence score for inclusion of a variant call from cortex
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
            max_snps_per_cluster=None, source='cluster_vcf_records', min_SNP_qual=None, min_dp4=None, min_GT_conf=None):
        self.vcf_files = vcf_files
        self.reference_seqs = {}
        pyfastaq.tasks.file_to_dict(reference_fasta, self.reference_seqs)
        for seq in self.reference_seqs.values():
            seq.seq = seq.seq.upper()
        self.vcf_outfile = vcf_outfile
        self.max_distance_between_variants = max_distance_between_variants
        self.homozygous_only = homozygous_only
        self.min_SNP_qual = min_SNP_qual
        self.min_dp4 = min_dp4
        self.min_GT_conf = min_GT_conf
        self.source = source
        self.max_REF_len = max_REF_len
        self.max_snps_per_cluster = max_snps_per_cluster


    @classmethod
    def _load_vcf_files(cls, filename_list, homozygous_only=False, max_REF_len=None, min_SNP_qual=None, min_dp4=None, min_GT_conf=None):
        '''Loads all the vcf files from filename_list. Returns tuple of:
        1. Sample name. If more than one sample name found, uses the first one
        and warns to stderr
        2. Dictionary. filename => list of header lines for that file
        3. Dictionary. ref name => list of VcfRecords sorted by position'''
        headers = {}
        vcf_records = None
        sample_name = None

        for filename in filename_list:
            headers[filename], new_records = vcf_file_read.vcf_file_to_dict(filename, homozygous_only=homozygous_only, remove_asterisk_alts=True, max_REF_len=max_REF_len, remove_useless_start_nucleotides=True, min_SNP_qual=min_SNP_qual, min_dp4=min_dp4, min_GT_conf=min_GT_conf)

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
    def _cluster_vcf_record_list(cls, vcf_records, max_distance_between_variants=1):
        
        indel_positions = []            #intialize variables and assign indel positions and lengths
        indel_lengths = {}
        for vcf_record in vcf_records:
            if vcf_record.INFO.get('SVTYPE') == 'INS':
                indel_positions.append(vcf_record.POS)
                indel_lengths[vcf_record.POS] = len(vcf_record.ALT)
            elif vcf_record.INFO.get('SVTYPE') == 'DEL':
                indel_postitions.append(vcf_record.POS)
                indel_lengths[vcf_record.POS] = len(vcf_record.ALT)
            else:
                pass

        count=0                     #delete vcf records that lie within indel regions
        for vcf_record in vcf_records:
            if len(indel_positions) > 1:
                for i in len(indel_positions):
                    if vcf_record.POS >= indel_position[i] and vcf_record.POS <= (indel_position[i]+indel_lengths[i]) and vcf_record.QUAL != None:
                        vcf_records.pop(count)
                    elif vcf_record.POS >= indel_position[i] and vcf_record.POS <= (indel_position[i]+indel_lengths[i]) and vcf_record.QUAL == None:
                        print('\nWarning: Tossing out overlapping cortex indel calls.\n')
                        vcf_records.pop(count)
                    else:
                        pass
            elif vcf_record.POS >= indel_positions[0] and vcf_record.POS <= (indel_positions[0]+indel_lengths[indel_positions[0]]) and vcf_record.QUAL != None:
                vcf_records.pop(count)
            elif vcf_record.POS >= indel_positions[0] and vcf_record.POS <= (indel_positions[0]+indel_lengths[indel_positions[0]]) and vcf_record.QUAL == None:
                print('\nWarning: Tossing out overlapping cortex indel calls.\n')
                vcf_records.pop(count)
            else:
                pass
            count+=1
                        
        new_list = [vcf_record_cluster.VcfRecordCluster(max_distance_between_variants=max_distance_between_variants)]
        
        for vcf_record in vcf_records:
            if not new_list[-1].add_vcf_record(vcf_record):
                            new_list.append(vcf_record_cluster.VcfRecordCluster(vcf_record=vcf_record, max_distance_between_variants=max_distance_between_variants))
        

                            
        return new_list


    def run(self):
        sample_name, vcf_headers, vcf_records = VcfClusterer._load_vcf_files(self.vcf_files, homozygous_only=self.homozygous_only, max_REF_len=self.max_REF_len, min_SNP_qual=self.min_SNP_qual, min_dp4=self.min_dp4, min_GT_conf=self.min_GT_conf)

        f_out = pyfastaq.utils.open_file_write(self.vcf_outfile)
        print('##fileformat=VCFv4.2', file=f_out)
        print('##source=', self.source, sep='', file=f_out)
        print('#CHROM', 'POS', 'ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name, sep='\t', file=f_out)

        for ref_name in sorted(vcf_records):
            ref_seq = self.reference_seqs[ref_name]
            cluster_list = VcfClusterer._cluster_vcf_record_list(vcf_records[ref_name], max_distance_between_variants=self.max_distance_between_variants)
            
            for cluster in cluster_list:
                clustered_vcf = cluster.make_one_merged_vcf_record_for_gramtools(ref_seq, max_snps=self.max_snps_per_cluster)
                if clustered_vcf is not None:
                    print(clustered_vcf, file=f_out)

        pyfastaq.utils.close(f_out)
