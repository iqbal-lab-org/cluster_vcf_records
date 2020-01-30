import typing
import logging
import operator

import pyfastaq

from cluster_vcf_records import vcf_file_read, vcf_record_cluster


class VcfClusterer:
    """Class to cluster one (or more) VCF files. Records in the VCf files
    can be in any order.

    Required parameters:
        vcf_files: list of VCF files
        reference_fasta: FASTA file of reference genome. Must be the same one
           that was used to make the VCF files
        vcf_outfile: name of output VCF file

    Optional parameters:
        cluster_boundary_size: Any variants that are within this distance of a cluster
           start & end positional boundaries will be put in that cluster.
           The default of 0 means that only overlapping variants are clustered.
        homozygous_only: Set this to True to only load homozygous variants from
           the input VCF files, ie where the genotype is 1/1.
        max_REF_len: When loading the VCF files, any records with REF longer
           than this will be ignored. By default, there is no maximum.
           This option should not normally be needed. It was intended to
           liimit the number of alleles, so use max_alleles_per_cluster instead
        max_alleles_per_cluster: maximum allowed alleles in one cluster. If a cluster
           has more than this, then combintations of SNPs are not generated.
           Instead, each indel is used on its own to make ALTs,
           and all SNPs are applied to make another ALT with all the SNPs.
        source: this is put into the source=foo part of the header of the
           output VCF file.

    Use it like this:

    clusterer = VcfCluster(vcf_files, reference_fasta, vcf_outfile[, options...])
    clusterer.run()
        """

    def __init__(
        self,
        vcf_files,
        reference_fasta,
        vcf_outfile,
        cluster_boundary_size=0,
        homozygous_only=False,
        max_REF_len=None,
        max_alleles_per_cluster=None,
        source="cluster_vcf_records",
        merge_method="gramtools",
        max_gap_indel_rmdup=100,
    ):
        self.vcf_files = vcf_files
        self.reference_seqs = {}
        pyfastaq.tasks.file_to_dict(reference_fasta, self.reference_seqs)

        # If fasta header contains whitespace, add the first word as valid sequence ID (on top of whole line- pyfastaq module produced).
        # This allows vcf files with CHROM field being only first word of fasta header, to be parsed correctly.
        stripped_keys = {
            key.split()[0]: self.reference_seqs[key] for key in self.reference_seqs
        }
        self.reference_seqs.update(stripped_keys)

        for seq in self.reference_seqs.values():
            seq.seq = seq.seq.upper()
        self.vcf_outfile = vcf_outfile
        self.cluster_boundary_size = cluster_boundary_size
        self.homozygous_only = homozygous_only
        self.source = source
        self.max_REF_len = max_REF_len
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.merge_method = merge_method
        self.max_gap_indel_rmdup = max_gap_indel_rmdup

        allowed_merge_methods = {"gramtools", "simple", "gt_aware"}
        if self.merge_method not in {"gramtools", "simple", "gt_aware"}:
            raise RuntimeError(
                'Error! merge_method "'
                + self.merge_method
                + '" not allowed. Must be one of: '
                + ",".join(sorted(list(allowed_merge_methods)))
            )

    @classmethod
    def _load_vcf_files(
        cls,
        filename_list,
        reference_seqs,
        homozygous_only=False,
        max_REF_len=None,
        min_SNP_qual=None,
        min_dp4=None,
        min_GT_conf=None,
    ):
        """Loads all the vcf files from filename_list. Returns tuple of:
        1. Sample name. If more than one sample name found, uses the first one
        and warns to stderr
        2. Dictionary. filename => list of header lines for that file
        3. Dictionary. ref name => list of VcfRecords sorted by position.

        reference_seqs should be a dictionary of sequence name -> sequence.
        This causes all records from the VCF to be sanity checked against the reference sequence,
        and any records where the REF seq does not match the expected sequence is removed."""
        headers = {}
        vcf_records = None
        sample_name = None

        for filename in filename_list:
            headers[filename], new_records = vcf_file_read.vcf_file_to_dict(
                filename,
                homozygous_only=homozygous_only,
                remove_asterisk_alts=True,
                max_REF_len=max_REF_len,
                remove_useless_start_nucleotides=True,
                min_SNP_qual=min_SNP_qual,
                min_dp4=min_dp4,
                min_GT_conf=min_GT_conf,
                reference_seqs=reference_seqs,
                error_on_bad_POS=False,
            )

            new_sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(
                headers[filename]
            )
            if sample_name is None and new_sample_name is not None:
                sample_name = new_sample_name
            elif new_sample_name != sample_name:
                logging.warning(
                    'Using first sample name found "'
                    + str(sample_name)
                    + '". Found a different (or no) sample name "'
                    + str(new_sample_name)
                    + '", which will not be used'
                )

            if vcf_records is None:
                vcf_records = new_records
            else:
                for ref_name, record_list in new_records.items():
                    if ref_name not in vcf_records:
                        vcf_records[ref_name] = record_list
                    else:
                        vcf_records[ref_name].extend(record_list)

        for record_list in vcf_records.values():
            record_list.sort(key=operator.attrgetter("POS"))

        if sample_name is None:
            logging.warning('No sample name found in VCF files. Going to use "sample"')
            sample_name = "sample"

        return sample_name, headers, vcf_records

    @classmethod
    def _expand_alts_in_vcf_record_list(cls, vcf_records):
        """Input: list of vcf_records. Returns new list, where
        any records with >ALT is replaced with one vcf record per ALT.
        This doesn't change FORMAT or INFO columns, which means they
        are now broken for those records"""
        new_vcf_records = []
        for record in vcf_records:
            new_vcf_records.extend(record.to_record_per_alt())
        return new_vcf_records

    @classmethod
    def _expand_alts_and_remove_duplicates_in_list(
        cls, vcf_records, ref_seq, indel_gap=100
    ):
        """Input: list of VCF records, all from the same CHROM. ref_seq = sequence
        of that CHROM. Expands any record in the list that has >ALT, into
        one record per ALT. Removes duplicated records, where REF and ALT
        are the same (at the same position!), or where there is the same
        indel more than once, but written in a different way (eg indel in
        homopolymer run can be put in >1 way in a VCF. Checks indels
        are the same within indel_gap nucleotides of each other"""
        expanded_vcf_records = VcfClusterer._expand_alts_in_vcf_record_list(vcf_records)
        new_vcf_records = [x for x in expanded_vcf_records if not x.is_snp()]

        for i in range(len(new_vcf_records) - 1):
            j = i + 1
            while (
                j < len(new_vcf_records)
                and new_vcf_records[i].ref_end_pos() + indel_gap
                > new_vcf_records[j].POS
            ):
                if new_vcf_records[i].is_the_same_indel(new_vcf_records[j], ref_seq):
                    new_vcf_records.pop(j)
                else:
                    j += 1

        new_vcf_records.extend([x for x in expanded_vcf_records if x.is_snp()])
        new_vcf_records.sort(key=operator.attrgetter("POS"))
        return new_vcf_records

    @classmethod
    def _cluster_vcf_record_list(cls, vcf_records, cluster_boundary_size=0):
        new_cluster_list = [
            vcf_record_cluster.VcfRecordCluster(
                cluster_boundary_size=cluster_boundary_size
            )
        ]

        # We try adding each vcf_record to the lastmost cluster; if this fails, we put it in a new cluster of its own.
        for vcf_record in vcf_records:
            last_cluster = new_cluster_list[-1]
            successfully_added = last_cluster.add_vcf_record(vcf_record)

            if not successfully_added:  # Make a new cluster
                new_cluster = vcf_record_cluster.VcfRecordCluster(
                    vcf_record=vcf_record, cluster_boundary_size=cluster_boundary_size,
                )
                new_cluster_list.append(new_cluster)

        return new_cluster_list

    def run(self):
        sample_name, vcf_headers, vcf_records = VcfClusterer._load_vcf_files(
            self.vcf_files,
            self.reference_seqs,
            homozygous_only=self.homozygous_only,
            max_REF_len=self.max_REF_len,
        )

        f_out = pyfastaq.utils.open_file_write(self.vcf_outfile)
        print("##fileformat=VCFv4.2", file=f_out)
        print("##source=", self.source, sep="", file=f_out)
        print(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            sample_name,
            sep="\t",
            file=f_out,
        )

        for ref_name in sorted(vcf_records):
            ref_seq = self.reference_seqs[ref_name]

            if self.merge_method == "gramtools":
                rmdup_list = VcfClusterer._expand_alts_and_remove_duplicates_in_list(
                    vcf_records[ref_name], ref_seq, indel_gap=self.max_gap_indel_rmdup
                )

                cluster_list = VcfClusterer._cluster_vcf_record_list(
                    rmdup_list, cluster_boundary_size=self.cluster_boundary_size,
                )

                for cluster in cluster_list:
                    if len(cluster) > 0:
                        clustered_vcf = cluster.make_one_merged_vcf_record_for_gramtools(
                            ref_seq, max_alleles=self.max_alleles_per_cluster
                        )
                        if clustered_vcf is not None:
                            print(clustered_vcf, file=f_out)
                        else:
                            merged_record = cluster.make_separate_indels_and_one_alt_with_all_snps_no_combinations(
                                ref_seq
                            )
                            if merged_record is not None:
                                print(merged_record, file=f_out)

            elif self.merge_method == "simple":
                cluster_list = VcfClusterer._cluster_vcf_record_list(
                    vcf_records[ref_name],
                    cluster_boundary_size=self.cluster_boundary_size,
                )
                for cluster in cluster_list:
                    clustered_vcf = cluster.make_simple_merged_vcf_with_no_combinations(
                        ref_seq
                    )
                    for vcf in cluster.vcf_records:
                        print(vcf, file=f_out)
            elif self.merge_method == "gt_aware":
                cluster_list = VcfClusterer._cluster_vcf_record_list(
                    vcf_records[ref_name],
                    cluster_boundary_size=self.cluster_boundary_size,
                )
                for cluster in cluster_list:
                    clustered_vcf = cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(
                        ref_seq
                    )
                    for vcf in cluster.vcf_records:
                        print(vcf, file=f_out)
            else:
                raise RuntimeError(
                    'merge_method "'
                    + self.merge_method
                    + '" not recognised. Cannot continue'
                )

        pyfastaq.utils.close(f_out)


def cluster(
    input_vcf_file_paths: typing.List[str],
    reference_file_path: str,
    output_vcf_file_path: str,
    **kw
):
    _vcf_cluster = VcfClusterer(
        input_vcf_file_paths, reference_file_path, output_vcf_file_path, **kw
    )
    return _vcf_cluster.run()
