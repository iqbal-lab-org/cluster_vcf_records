import copy
import itertools
import logging

from cluster_vcf_records import vcf_record

class Error (Exception): pass

class VcfRecordCluster:
    def __init__(self, vcf_record=None, max_distance_between_variants=31):
        self.vcf_records = [] if vcf_record is None else [vcf_record]
        self.max_distance_between_variants = max_distance_between_variants


    def __getitem__(self, i):
        return self.vcf_records[i]


    def __len__(self):
        return len(self.vcf_records)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def add_vcf_record(self, vcf_record):
        if len(self) == 0 or True in {x.near_to_position(vcf_record.POS, self.max_distance_between_variants) for x in self.vcf_records}:
            self.vcf_records.append(vcf_record)
            return True
        else:
            return False


    def start_and_end(self):
        if len(self) == 0:
            return None, None
        else:
            return min([x.POS for x in self.vcf_records]), max([x.ref_end_pos() for x in self.vcf_records])


    def make_one_merged_vcf_record_for_gramtools(self, ref_seq, max_snps=None):
        '''Returns one new VcfRecord that can be used as input to gramtools.
        It pads the reference if necessary, and lists all the variants
        (including all combinations of SNPs) in the ALT field of the
        VcfRecord.
        Note: gramtools needs PASS in the filter column, so the returned
        VcfRecord always has PASS.'''
        if len(self) == 0:
            return None
        elif len(self) == 1:
            record = copy.copy(self[0])
            record.FILTER = 'PASS'
            return record

        logging.debug('make_one_merged_vcf_record_for_gramtools() start. Number of records: ' +  str(len(self)))
        for record in self.vcf_records:
            logging.debug('make_one_merged_vcf_record_for_gramtools() input record: ' + str(record))

        # Gather together the SNP and non-SNP alleles.
        # Also sanity check that the CHROM names are all the same
        # and determine the final start and end positions of the
        # vcf record we will output
        nucleotides = {'A', 'C', 'G', 'T'}
        snps = {} # position => set of alts (and the ref nucleotide)
        non_snps = [] # list of tuples (ref position, ref seq, alt seq)
        chrom_names = set()
        final_start = float('Inf')
        final_end = -1

        for record in self.vcf_records:
            final_start = min(final_start, record.POS)
            final_end = max(final_end, record.ref_end_pos())
            chrom_names.add(record.CHROM)

            if record.REF in nucleotides:
                for alt in record.ALT:
                    if alt in nucleotides:
                        if record.POS not in snps:
                            snps[record.POS] = {record.REF}
                        snps[record.POS].add(alt)
                    else:
                        non_snps.append((record.POS, record.REF, alt))
            else:
                for alt in record.ALT:
                    non_snps.append((record.POS, record.REF, alt))

        if len(chrom_names) != 1:
            raise Error('Error! More than one CHROM found. Got:' + str(chrom_names))
        chrom_name = chrom_names.pop()

        # generate all the allele combinations from the SNPs.
        snp_positions = []
        snp_nucleotides = []
        for position in sorted(snps):
            snp_positions.append(position)
            snp_nucleotides.append(sorted(list(snps[position])))
        ref_seq_for_vcf = ref_seq[final_start:final_end+1]
        alleles = set()

        logging.debug('make_one_merged_vcf_record_for_gramtools() alleles from SNPs: ' + str(snp_nucleotides))
        if max_snps is not None and len(snp_nucleotides) > max_snps:
            logging.info('Skip cluster because too many (' + str(len(snp_nucleotides)) + ') SNPs. Max allowed is ' + str(max_snps))
            for record in self.vcf_records:
                logging.info('    SKIP RECORD: ' + str(record))
            return None

        for combination in itertools.product(*snp_nucleotides):
            alt_seq = list(ref_seq_for_vcf)
            for i, position in enumerate(snp_positions):
                alt_seq[position - final_start] = combination[i]
            alleles.add(''.join(alt_seq))

        # remove the ref allele, because it goes in the REF
        # column of the VCF
        alleles.remove(ref_seq_for_vcf)

        logging.debug('make_one_merged_vcf_record_for_gramtools() generate non-SNP alleles')
        # add in the non-snp alleles. Need to add the flanking regions
        for pos, ref, alt in non_snps:
            fields = [chrom_name, str(pos + 1), '.', ref, alt, '.', 'PASS', 'SVTYPE=COMPLEX']
            new_vcf = vcf_record.VcfRecord('\t'.join(fields))
            new_vcf.add_flanking_seqs(ref_seq, final_start, final_end)
            alleles.add(new_vcf.ALT[0])

        alleles = sorted(list(alleles))

        fields = [chrom_name, str(final_start + 1), '.', ref_seq_for_vcf,
                  ','.join(alleles), '.', 'PASS', 'SVTYPE=COMPLEX']
        logging.debug('make_one_merged_vcf_record_for_gramtools number of alts: ' + str(len(alleles)))
        return vcf_record.VcfRecord('\t'.join(fields))

