import copy
import itertools
import logging

from cluster_vcf_records import vcf_record


class VcfRecordCluster:
    def __init__(self, vcf_record=None, cluster_boundary_size=0):
        if vcf_record is None:
            self.vcf_records = []
            self.start = None
            self.end = None
        else:
            self.vcf_records = [vcf_record]
            self.start = vcf_record.POS
            self.end = self.start + len(vcf_record.REF) - 1

        self.cluster_boundary_size = cluster_boundary_size

    def __getitem__(self, i):
        return self.vcf_records[i]

    def __len__(self):
        return len(self.vcf_records)

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def add_vcf_record(self, vcf_record):
        """
        Query the current cluster boundaries allowing for a distance around them, and include vcf_record if it is inside.
        """
        assert vcf_record.POS >= 0  # constructor should enforce this anyway
        record_end = vcf_record.POS + len(vcf_record.REF) - 1
        if len(self) == 0:
            self.vcf_records.append(vcf_record)
            self.start = vcf_record.POS
            self.end = record_end
            return True

        else:
            # Here not making assumption that vcf_record's POS is >= self.start; yet for records processed in sorted order,
            # this should be the case.
            if (
                self.start - self.cluster_boundary_size
                <= vcf_record.POS
                <= self.end + self.cluster_boundary_size
            ):
                self.vcf_records.append(vcf_record)

                # Update the cluster boundaries if necessary
                if record_end > self.end:
                    self.end = record_end
                elif vcf_record.POS < self.start:
                    self.start = vcf_record.POS
                return True

        return False

    def start_and_end(self):
        if len(self) == 0:
            return None, None
        else:
            return (
                min([x.POS for x in self.vcf_records]),
                max([x.ref_end_pos() for x in self.vcf_records]),
            )

    def make_one_merged_vcf_record_for_gramtools(self, ref_seq, max_alleles=5000):
        """Returns one new VcfRecord that can be used as input to gramtools.
        It pads the reference if necessary, and lists all the variants
        (including all combinations of SNPs) in the ALT field of the
        VcfRecord.
        Note: gramtools needs PASS in the filter column, so the returned
        VcfRecord always has PASS."""
        if len(self) == 0:
            return None
        elif len(self) == 1:
            record = copy.copy(self[0])
            record.FILTER = {"PASS"}
            return record

        logging.debug(
            "make_one_merged_vcf_record_for_gramtools() start. Number of records: "
            + str(len(self))
        )
        for record in self.vcf_records:
            logging.debug(
                "make_one_merged_vcf_record_for_gramtools() input record: "
                + str(record)
            )

        # Gather together the SNP and non-SNP alleles.
        # Also sanity check that the CHROM names are all the same
        # and determine the final start and end positions of the
        # vcf record we will output
        nucleotides = {"A", "C", "G", "T"}
        snps = {}  # position => set of alts (and the ref nucleotide)
        non_snps = []  #  list of tuples (ref position, ref seq, alt seq)
        chrom_names = set()
        final_start = float("Inf")
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
            raise Exception("Error! More than one CHROM found. Got:" + str(chrom_names))
        chrom_name = chrom_names.pop()

        # generate all the allele combinations from the SNPs.
        snp_positions = []
        snp_nucleotides = []
        for position in sorted(snps):
            snp_positions.append(position)
            snp_nucleotides.append(sorted(list(snps[position])))
        ref_seq_for_vcf = ref_seq[final_start : final_end + 1]
        alleles = set()

        # work out min total alleles without making them. Making them could
        #  take a long time if too many! Can onnly put lower bound on
        # the final unique number because all the combinations may have duplicates.
        total_alleles_lower_bound = 1
        for x in snp_nucleotides:
            total_alleles_lower_bound *= len(x)
        total_alleles_lower_bound += len(non_snps)

        if max_alleles is not None and total_alleles_lower_bound > max_alleles:
            return None

        for combination in itertools.product(*snp_nucleotides):
            alt_seq = list(ref_seq_for_vcf)
            for i, position in enumerate(snp_positions):
                # In a VCF file, POS is 1-based.
                # We're storing POS in memory 0-based.
                # Therefore we should never have position < 0.
                # The VcfRecord constructor should prevent this, but check anyway
                if position < 0:
                    raise ValueError(
                        f"POS value ({position+1}) less than 1. Ensure VCF is valid."
                    )
                alt_seq[position - final_start] = combination[i]
            alleles.add("".join(alt_seq))
            for non_snp_pos, non_snp_ref, non_snp_alt in non_snps:
                start_pos = non_snp_pos - final_start
                new_seq = (
                    alt_seq[:start_pos]
                    + list(non_snp_alt)
                    + alt_seq[start_pos + len(non_snp_ref) :]
                )
                alleles.add("".join(new_seq))

        # remove the ref allele (if it's there), because it goes in the REF
        # column of the VCF
        try:
            alleles.remove(ref_seq_for_vcf)
        except:
            pass

        if max_alleles is not None and len(alleles) > max_alleles:
            return None
        alleles = sorted(list(alleles))
        fields = [
            chrom_name,
            str(final_start + 1),
            ".",
            ref_seq_for_vcf,
            ",".join(alleles),
            ".",
            "PASS",
            "SVTYPE=COMPLEX",
        ]
        logging.debug(
            "make_one_merged_vcf_record_for_gramtools number of alts: "
            + str(len(alleles))
        )
        return vcf_record.VcfRecord("\t".join(fields))

    def make_simple_merged_vcf_with_no_combinations(self, ref_seq):
        """Does a simple merging of all variants in this cluster.
        Assumes one ALT in each variant. Uses the ALT for each
        variant, making one new vcf_record that has all the variants
        put together"""
        if len(self) <= 1:
            return

        merged_vcf_record = self.vcf_records[0]

        for i in range(1, len(self.vcf_records), 1):
            if self.vcf_records[i].intersects(merged_vcf_record):
                return
            else:
                merged_vcf_record = merged_vcf_record.merge(
                    self.vcf_records[i], ref_seq
                )

        self.vcf_records = [merged_vcf_record]

    def make_simple_gt_aware_merged_vcf_with_no_combinations(self, ref_seq):
        """Does a simple merging of all variants in this cluster.
        Assumes one ALT in each variant. Uses the called allele for each
        variant, making one new vcf_record that has all the variants
        put together"""
        if len(self) <= 1:
            return

        merged_vcf_record = self.vcf_records[0]

        for i in range(1, len(self.vcf_records), 1):
            if self.vcf_records[i].intersects(merged_vcf_record):
                return
            else:
                merged_vcf_record = merged_vcf_record.gt_aware_merge(
                    self.vcf_records[i], ref_seq
                )

        self.vcf_records = [merged_vcf_record]

    def make_separate_indels_and_one_alt_with_all_snps_no_combinations(self, ref_seq):
        """Returns a VCF record, where each indel from this
        cluster is in a separate ALT. Then all the remaining SNPs are
        applied to make one ALT. If >1 SNP in same place, either one
        might be used"""
        final_start_position = min([x.POS for x in self.vcf_records])
        final_end_position = max([x.ref_end_pos() for x in self.vcf_records])
        snps = []
        new_vcf_records = []

        for record in self.vcf_records:
            if record.is_snp():
                snps.append(copy.copy(record))
            else:
                new_record = copy.copy(record)
                new_record.add_flanking_seqs(
                    ref_seq, final_start_position, final_end_position
                )
                new_vcf_records.append(new_record)

        if len(snps):
            new_record = copy.copy(snps[0])
            for snp in snps[1:]:
                merged = new_record.merge(snp, ref_seq)
                if merged is not None:
                    new_record = merged
            new_record.add_flanking_seqs(
                ref_seq, final_start_position, final_end_position
            )
            new_vcf_records.append(new_record)

        alts = ",".join(sorted(list(set([x.ALT[0] for x in new_vcf_records]))))
        new_record = vcf_record.VcfRecord(
            "\t".join(
                [
                    self.vcf_records[0].CHROM,
                    str(final_start_position + 1),
                    ".",
                    new_vcf_records[0].REF,
                    alts,
                    ".",
                    "PASS",
                    ".",
                ]
            )
        )

        return new_record
