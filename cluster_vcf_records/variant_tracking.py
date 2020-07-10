from collections import namedtuple
import itertools
import json
import logging
import multiprocessing
from operator import attrgetter
import os
import re
import tempfile

from bitarray import bitarray
import pyfastaq
import pysam

from cluster_vcf_records import allele_combinations, utils, vcf_file_read, vcf_record


Variant = namedtuple("Variant", ["seq_id", "pos", "ref", "alt"])


def _run_one_cluster_chunk(tracker, outprefix, max_ref_length, max_alleles, index_range):
    logging.info(f"Start clustering chunk {index_range[0]}-{index_range[1]}")
    tracker.cluster_one_chunk(outprefix, max_ref_length, max_alleles=max_alleles, split_start=index_range[0], split_end=index_range[1])



def vcf_records_make_same_allele_combination(record1, record2, ref_seqs):
    if record1.CHROM != record2.CHROM:
        return False
    assert record1.ref_end_pos() < record2.POS
    alleles1 = record1.inferred_var_seqs_plus_flanks(ref_seqs[record1.CHROM], 0)[1]
    alleles2 = record2.inferred_var_seqs_plus_flanks(ref_seqs[record2.CHROM], 0)[1]
    join_seq = ref_seqs[record1.CHROM][record1.ref_end_pos() + 1 : record2.POS]

    combinations = set()
    for a1 in alleles1:
        for a2 in alleles2:
            allele = a1 + join_seq + a2
            if allele in combinations:
                return True
            combinations.add(allele)
    return False


def variants_overlap(var1, var2):
    end1 = var1.pos + len(var1.ref) - 1
    end2 = var2.pos + len(var2.ref) - 1
    return var1.seq_id == var2.seq_id and var1.pos <= end2 and var2.pos <= end1


def _load_one_vcf_file(
    vcf_file, ref_seqs, ref_seq_to_id, ref_fasta, temp_dir, break_alleles
):
    sample = vcf_file_read.get_sample_name_from_vcf_file(vcf_file)
    if sample is None:
        raise Exception(f"Error getting sample name from vcf file {vcf_file}")
    tmpdir = tempfile.mkdtemp(prefix="normalize_vcf.", dir=temp_dir)
    simplified_vcf = os.path.join(tmpdir, "simplified.vcf")
    normalized_vcf = os.path.join(tmpdir, "normalized.vcf")
    utils.simplify_vcf(vcf_file, simplified_vcf, ref_seqs=ref_seqs)
    utils.normalise_vcf(
        simplified_vcf, ref_fasta, normalized_vcf, break_alleles=break_alleles
    )
    os.unlink(simplified_vcf)
    variants = []

    with open(normalized_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            record = vcf_record.VcfRecord(line)
            if record.POS < 0:
                logging.warning(
                    f"VCF record with negative POS in file {vcf_file}. Ignoring: {record}"
                )
                continue
            elif record.CHROM not in ref_seqs:
                logging.warning(
                    f"CHROM not recognised in VCF record in file {vcf_file}. Ignoring: {record}"
                )
                continue
            elif not record.ref_string_matches_ref_sequence(ref_seqs[record.CHROM]):
                logging.warning(
                    f"REF string does not match reference seq in file {vcf_file}. Ignoring: {record}"
                )
                continue
            elif "GT" not in record.FORMAT:
                logging.warning(
                    f"No GT in VCF record in file {vcf_file}. Ignoring: {record}"
                )
                continue

            gt_indexes = re.split("[/|]", record.FORMAT["GT"])
            if "." in gt_indexes:
                continue
            gt_indexes = set([int(x) for x in gt_indexes])
            if gt_indexes == {0}:
                continue

            for i in gt_indexes:
                if i > 0:
                    variants.append(
                        Variant(
                            ref_seq_to_id[record.CHROM],
                            record.POS,
                            record.REF,
                            record.ALT[i - 1],
                        )
                    )

    utils.rm_rf(tmpdir)
    return sample, variants


class Variants:
    def __init__(self):
        self.vars = {}

    def __len__(self):
        return len(self.vars)

    def __getitem__(self, x):
        return self.vars[x]

    def __contains__(self, x):
        return x in self.vars

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def add(self, new_variant):
        if new_variant not in self.vars:
            self.vars[new_variant] = len(self.vars)
        return self.vars[new_variant]

    def sorted_iter(self):
        for var in sorted(self.vars.keys(), key=attrgetter("seq_id", "pos")):
            yield var, self.vars[var]

    def save_to_file(self, filename):
        with os.popen("bgzip -c > " + filename, "w") as f:
            for v, v_id in self.sorted_iter():
                print(v_id, v.seq_id, v.pos + 1, v.ref, v.alt, sep="\t", file=f)

    def load_from_file(self, filename):
        with os.popen("gunzip -c " + filename) as f:
            for line in f:
                var_id, seq_id, pos, ref, alt = line.rstrip().split("\t")
                var = Variant(int(seq_id), int(pos) - 1, ref, alt)
                self.vars[var] = int(var_id)


class VariantBlock:
    def __init__(self, number_of_variants=0):
        self.bitarrays = []
        self.add_variants(number_of_variants)

    def clear_samples(self):
        for array in self.bitarrays:
            del array[:]

    def number_of_samples(self):
        if len(self.bitarrays) == 0:
            return 0
        else:
            return self.bitarrays[0].length()

    def number_of_variants(self):
        return len(self.bitarrays)

    def add_variants(self, number_of_variants):
        b = self.number_of_samples() * bitarray("0")
        for i in range(number_of_variants):
            self.bitarrays.append(bitarray(b))

    def add_samples(self, number_of_samples):
        assert len(self.bitarrays) > 0
        for array in self.bitarrays:
            array.extend(number_of_samples * bitarray("0"))

    def set_variant(self, variant_index, sample_index):
        self.bitarrays[variant_index][sample_index] = True

    def has_variant(self, variant_index, sample_index):
        return self.bitarrays[variant_index][sample_index]

    def write_to_bgzip_file_and_tab_index(self, outfile, variants):
        """Writes a bgzipped and tabix indexed file. The pos column is
        1-based and does not include end points"""
        logging.info(f"Saving batch of data to disk {outfile}")
        with os.popen("bgzip -c > " + outfile, "w") as f:
            print("#seq_id\tpos\tvar_id\tarray", file=f)
            for var, var_id in variants.sorted_iter():
                # No point writing the array if none of the samples have the variant
                if self.bitarrays[var_id].any():
                    print(
                        var.seq_id,
                        var.pos + 1,
                        var_id,
                        self.bitarrays[var_id].to01(),
                        sep="\t",
                        file=f,
                    )

        pysam.tabix_index(outfile, seq_col=0, start_col=1, end_col=1, zerobased=True)

    def estimate_mem_in_gb_if_more_added(self, samples, variants):
        return (
            (samples + self.number_of_samples())
            * (variants + self.number_of_variants())
            / 8e9
        )


def load_slice_of_block(infile, wanted_variants, seq_id, start, end):
    """Returns dictionary of variant id -> bitarray.
    start and end should be zero-based and include end points.
    wanted_variants should contain the variant ids of the variants to return.
    Any variant found in the block that is not in wanted_variants is ignored."""
    variants = {}
    if not isinstance(infile, pysam.libctabix.TabixFile):
        infile = pysam.TabixFile(infile)

    for line in infile.fetch(str(seq_id), start, end + 1):
        _, _, var_id, array = line.rstrip().split()
        var_id = int(var_id)
        if var_id in wanted_variants:
            variants[int(var_id)] = bitarray(array)
    return variants


def var_patterns_from_block_slices(block_files, wanted_variants, seq_id, start, end):
    """Returns a set of tuples. Each tuple is a combination of variants that is
    seen in one or more samples. Onlt the variants with ids in wanted_variants
    are used. All other variants are ignored"""
    var_patterns = set()
    for tabix_file in block_files:
        try:
            block = load_slice_of_block(tabix_file, wanted_variants, seq_id, start, end)
        except ValueError:  # happens if region is not in tabix file
            continue
        if len(block) == 0:
            continue
        sorted_vars = sorted(list(block.keys()))
        samples = block[sorted_vars[0]].length()

        for i in range(samples):
            var_pattern = tuple(sorted([v for v in sorted_vars if block[v][i]]))
            if len(var_pattern) > 0:
                var_patterns.add(var_pattern)

    return var_patterns


def var_pattern_to_alleles(
    variants, var_pattern, ref_seq, start, end, max_alleles=None
):
    sorted_vars = sorted(
        [variants[v] for v in var_pattern], key=attrgetter("seq_id", "pos")
    )
    allele = list(ref_seq[start : end + 1])
    for v1, v2 in itertools.combinations(sorted_vars, 2):
        if variants_overlap(v1, v2):
            (
                this_start,
                this_end,
                alts,
            ) = allele_combinations.var_cluster_to_coords_and_alts(
                sorted_vars, ref_seq, max_alleles=max_alleles
            )
            if alts is not None and start <= this_start <= this_end <= end:
                alts = {
                    ref_seq[start:this_start] + x + ref_seq[this_end + 1 : end + 1]
                    for x in alts
                }
                return alts
            else:
                return None

    # Do in reverse order so indels don't mess up coords
    for var in reversed(sorted_vars):
        assert start <= var.pos <= var.pos + len(var.ref) - 1 <= end
        allele_pos = var.pos - start
        allele[allele_pos : allele_pos + len(var.ref)] = var.alt

    return {"".join(allele)}


class VariantTracker:
    def __init__(self, root_dir, ref_fasta, mem_limit=2):
        self.root_dir = os.path.abspath(root_dir)
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.mem_lmit = mem_limit
        self.variants = Variants()
        self.cluster_limit = 8
        self.max_var_gap = 10
        (
            self.ref_seqs,
            self.ref_seq_names,
            self.ref_seq_to_id,
        ) = VariantTracker.load_ref_seq_data(self.ref_fasta)
        self.metadata_file = os.path.join(root_dir, "metadata.json")
        self.variants_file = os.path.join(root_dir, "variants.tsv.gz")
        if os.path.exists(root_dir):
            self._load_data_from_root()

    @classmethod
    def load_ref_seq_data(cls, fasta_file):
        seqs = {}
        pyfastaq.tasks.file_to_dict(fasta_file, seqs)
        ref_seq_names = sorted(seqs.keys())
        ref_seq_to_id = {x: i for i, x in enumerate(ref_seq_names)}
        return seqs, ref_seq_names, ref_seq_to_id

    def _load_data_from_root(self):
        with open(self.metadata_file) as f:
            self.var_block_files, self.samples = json.load(f)
        self.variants.load_from_file(self.variants_file)
        self.var_block_tabixes = [
            os.path.join(self.root_dir, x) for x in self.var_block_files
        ]

    def _add_var_block_file(self):
        filename = f"block.{len(self.var_block_files)}.tsv.gz"
        self.var_block_files.append(filename)
        self.samples.append([])

    def _write_last_var_block_file(self, var_block):
        filename = os.path.join(self.root_dir, self.var_block_files[-1])
        if not os.path.exists(filename):
            var_block.write_to_bgzip_file_and_tab_index(filename, self.variants)

    def _save_metadata_file(self):
        with open(self.metadata_file, "w") as f:
            json.dump([self.var_block_files, self.samples], f, indent=2, sort_keys=True)

    def _block_will_break_limits(self, var_block, mem_limit, sample_limit):
        over_sample_limit = (
            sample_limit is not None and var_block.number_of_samples() >= sample_limit
        )
        over_mem_limit = var_block.estimate_mem_in_gb_if_more_added(1, 1) > mem_limit
        return var_block.number_of_samples() > 0 and (
            over_sample_limit or over_mem_limit
        )

    def merge_vcf_files(
        self,
        infiles,
        temp_dir=None,
        cpus=1,
        mem_limit=2,
        force=False,
        sample_limit=None,
        break_alleles=True,
    ):
        if force:
            utils.rm_rf(self.root_dir)
        os.mkdir(self.root_dir)

        if temp_dir is None:
            temp_dir = os.path.join(self.root_dir, "tmp")
        if os.path.exists(temp_dir):
            made_temp_dir = False
        else:
            os.mkdir(temp_dir)
            made_temp_dir = True

        self.var_block_files = []
        self.samples = []
        self._add_var_block_file()
        self.variants = Variants()
        var_block = VariantBlock()

        for i in range(0, len(infiles), cpus):
            with multiprocessing.Pool(cpus) as pool:
                new_variants_lists = (
                    self.var_block_filesew_variants_lists
                ) = pool.starmap(
                    _load_one_vcf_file,
                    zip(
                        infiles[i : i + cpus],
                        itertools.repeat(self.ref_seqs),
                        itertools.repeat(self.ref_seq_to_id),
                        itertools.repeat(self.ref_fasta),
                        itertools.repeat(temp_dir),
                        itertools.repeat(break_alleles),
                    ),
                )

            for sample, new_variants in new_variants_lists:
                if len(new_variants) == 0:
                    continue

                if self._block_will_break_limits(var_block, mem_limit, sample_limit):
                    self._write_last_var_block_file(var_block)
                    self._add_var_block_file()
                    var_block.clear_samples()

                self.samples[-1].append(sample)
                if var_block.number_of_samples() == 0 == var_block.number_of_variants():
                    var_block.add_variants(1)
                    self.variants.add(new_variants[0])

                var_block.add_samples(1)

                for variant in new_variants:
                    if variant not in self.variants:
                        var_block.add_variants(1)
                    var_id = self.variants.add(variant)
                    var_block.set_variant(var_id, -1)

            if i % 100 == 0:
                logging.info(f"Loaded {i+cpus} files out of {len(infiles)}")

        if made_temp_dir:
            utils.rm_rf(temp_dir)

        logging.info(f"Loaded all {len(infiles)} VCF files")
        self._write_last_var_block_file(var_block)

        logging.info(f"Saving metadata to file {self.metadata_file}")
        self._save_metadata_file()

        logging.info(f"Saving variants to file {self.variants_file}")
        self.variants.save_to_file(self.variants_file)
        logging.info("Finished merging VCF files")

    def _variant_cluster_to_vcf_line(self, variants, variant_ids, max_alleles=None):
        if len(variants) == 0:
            return None
        ref_seq = self.ref_seqs[self.ref_seq_names[variants[0].seq_id]]
        if logging.getLogger().level <= logging.DEBUG:
            logging.debug(f"Clustering variants:")
            for v in variants:
                logging.debug(f"  {ref_seq.id}\t{v.pos}\t{v.ref}\t{v.alt}")
        start, end, alts = allele_combinations.var_cluster_to_coords_and_alts(
            variants, ref_seq, max_alleles=max_alleles
        )
        info_field = "."
        if alts is None:
            alts = set()
            var_id_to_var = dict(zip(variant_ids, variants))
            var_patterns = var_patterns_from_block_slices(
                self.var_block_tabixes, var_id_to_var, variants[0].seq_id, start, end
            )
            for var_pattern in var_patterns:
                alt_alleles = var_pattern_to_alleles(
                    var_id_to_var,
                    var_pattern,
                    ref_seq,
                    start,
                    end,
                    max_alleles=max_alleles,
                )
                if alt_alleles is None:
                    logging.warning("Conflicting allele combination:")
                    for var_id in sorted(list(var_pattern)):
                        var = var_id_to_var[var_id]
                        logging.warning(
                            f"  {ref_seq.id} {var.pos+1} {var.ref} {var.alt}"
                        )
                else:
                    alts.update(alt_alleles)

            info_field = "High_variability"

        if len(alts) == 0:
            logging.warning("Could not make VCF record from these variants:")
            for variant in variants:
                logging.warning("    " + str(variant))
            return None
        else:
            return vcf_record.VcfRecord(
                "\t".join(
                    [
                        ref_seq.id,
                        str(start + 1),
                        ".",
                        ref_seq[start : end + 1],
                        ",".join(sorted(list(alts))),
                        ".",
                        "PASS",
                        info_field,
                    ]
                )
            )

    def _print_cluster_vcf_header(self, filehandle, max_alleles):
        print("##fileformat=VCFv4.2", file=filehandle)
        for seq in self.ref_seqs.values():
            print(f"##contig=<ID={seq.id},length={len(seq)}>", file=filehandle)
        print('##FILTER=<ID=PASS,Description="All filters passed">', file=filehandle)
        if max_alleles is not None:
            print(
                f'##INFO=<ID=High_variability,Number=0,Type=Flag,Description="Position is in a region of high variability, with more alleles than the limit of {max_alleles}">',
                file=filehandle,
            )
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=filehandle)

    def _init_clusters(self):
        self.clusters = []
        self.clusters_insertions = set()
        self.clusters_deletions = set()

    def _append_to_clusters(self, to_append):
        has_insertion = any([len(to_append[0].REF) < len(x) for x in to_append[0].ALT])
        has_deletion = any([len(to_append[0].REF) > len(x) for x in to_append[0].ALT])
        if has_insertion:
            self.clusters_insertions.add(len(self.clusters))
        if has_deletion:
            self.clusters_deletions.add(len(self.clusters))
        self.clusters.append(to_append)

    def _merge_last_records_in_clusters(self, max_alleles):
        if len(self.clusters) < 2 or len(self.clusters_insertions) == 0 or len(self.clusters_deletions) == 0:
            return

        best_merged = None
        best_merged_index = None
        merged_record, merged_vars, merged_var_ids = self.clusters[-1]
        for i in reversed(range(len(self.clusters) - 1)):
            if merged_record.CHROM != self.clusters[i][0].CHROM or self.clusters[i][0].ref_end_pos() + self.max_var_gap < merged_record.POS:
                break

            have_same_combos = vcf_records_make_same_allele_combination(self.clusters[i][0], merged_record, self.ref_seqs)
            if i == 0 and not have_same_combos:
                break

            merged_vars = self.clusters[i][1] + merged_vars
            merged_var_ids = self.clusters[i][2] + merged_var_ids
            merged_record = self._variant_cluster_to_vcf_line(
                merged_vars, merged_var_ids, max_alleles=max_alleles
            )
            if merged_record is None:
                break

            if have_same_combos:
                best_merged = (merged_record, merged_vars, merged_var_ids)
                best_merged_index = i

        if best_merged_index is not None:
            self.clusters_insertions = {x for x in self.clusters_insertions if x < best_merged_index}
            self.clusters_deletions = {x for x in self.clusters_deletions if x < best_merged_index}
            del self.clusters[best_merged_index:]
            self._append_to_clusters(best_merged)


    def _merge_last_records_in_clusters_old(self, max_alleles):
        while len(self.clusters) > 1:
            if vcf_records_make_same_allele_combination(
                self.clusters[-2][0], self.clusters[-1][0], self.ref_seqs
            ):
                new_variants = self.clusters[-2][1] + self.clusters[-1][1]
                new_variant_ids = self.clusters[-2][2] + self.clusters[-1][2]
                new_record = self._variant_cluster_to_vcf_line(
                    new_variants, new_variant_ids, max_alleles=max_alleles
                )
                if new_record is not None:
                    self.clusters.pop()
                    self.clusters[-1] = (new_record, new_variants, new_variant_ids)
                else:
                    break
            else:
                break

    def _print_first_clusters(self, filehandle):
        while len(self.clusters) > self.cluster_limit:
            print(self.clusters[0][0], file=filehandle)
            self.clusters.pop(0)
            self.clusters_insertions.discard(0)
            self.clusters_deletions.discard(0)
            self.clusters_insertions = {x-1 for x in self.clusters_insertions}
            self.clusters_deletions = {x-1 for x in self.clusters_deletions}


    def find_parallel_cluster_split_points(self, splits):
        targets = [x * int(len(self.variants) / splits) for x in range(1, splits)]
        good_count = 0
        previous_var = None
        split_ranges = []
        looking_for_split = False
        target_index = 0

        for i, (var, var_id) in enumerate(self.variants.sorted_iter()):
            if not looking_for_split and i > targets[target_index]:
                if len(var.ref) == len(var.alt) == 1:
                    looking_for_split = True
                    good_count = 1
                    previous_var = var
                    continue

            if looking_for_split:
                if len(var.ref) == len(var.alt) == 1:
                    if var.pos > previous_var.pos:
                        good_count += 1
                    previous_var = var
                else:
                    good_count = 0

                if good_count > self.cluster_limit + 2:
                    if len(split_ranges) == 0:
                        split_ranges.append((0, i))
                    else:
                        split_ranges.append((split_ranges[-1][-1]+1, i))
                    target_index += 1
                    if target_index >= len(targets):
                        break
                    looking_for_split = False

        split_ranges.append((split_ranges[-1][-1]+1, len(self.variants) - 1))
        return split_ranges


    def cluster(self, outprefix, max_ref_length, max_alleles=None, cpus=1):
        if cpus > 1:
            split_ranges = self.find_parallel_cluster_split_points(cpus)
            outprefixes = [f"{outprefix}.{i}" for i in range(len(split_ranges))]
            logging.info(f"Running clustering in parallel using {cpus} cpus")
            with multiprocessing.Pool(cpus) as pool:
                pool.starmap(
                    _run_one_cluster_chunk,
                    zip(
                        itertools.repeat(self),
                        outprefixes,
                        itertools.repeat(max_ref_length),
                        itertools.repeat(max_alleles),
                        split_ranges
                    )
                )
            utils.cat_vcfs([x + ".vcf" for x in outprefixes], f"{outprefix}.vcf", delete_files=True)
            utils.cat_tsvs([x + ".excluded.tsv" for x in outprefixes], f"{outprefix}.excluded.tsv",  delete_files=True)
        else:
            self.cluster_one_chunk(outprefix, max_ref_length, max_alleles=max_alleles)


    def cluster_one_chunk(self, outprefix, max_ref_length, max_alleles=None, split_start=None, split_end=None):
        out_main = f"{outprefix}.vcf"
        out_exclude = f"{outprefix}.excluded.tsv"
        total_vars = len(self.variants) if split_start is None else split_end - split_start + 1

        with open(out_main, "w") as f_main, open(out_exclude, "w") as f_exclude:
            self._print_cluster_vcf_header(f_main, max_alleles)
            print("CHROM\tPOS\tREF\tALT", file=f_exclude)
            variants = []
            variant_ids = []
            var_count = 0
            report_count = int(total_vars / 10)
            self._init_clusters()

            for i, (var, var_id) in enumerate(self.variants.sorted_iter()):
                if split_start is not None:
                    if i < split_start:
                        continue
                    elif i > split_end:
                        break

                var_count += 1
                if var_count >= report_count:
                    pc_complete = round(100 * var_count / total_vars)
                    if split_start is None:
                        range_string = ""
                    else:
                        range_string = f" for split range {split_start}-{split_end}"
                    logging.info(
                        f"Completed {pc_complete}% ({var_count}/{total_vars}) of variants{range_string}"
                    )
                    report_count += int(total_vars / 10)

                if len(var.ref) > max_ref_length:
                    print(
                        self.ref_seq_names[var.seq_id],
                        var.pos + 1,
                        var.ref,
                        var.alt,
                        sep="\t",
                        file=f_exclude,
                    )
                    continue

                if len(variants) == 0 or any(
                    [variants_overlap(var, x) for x in variants]
                ):
                    variants.append(var)
                    variant_ids.append(var_id)
                    continue

                record = self._variant_cluster_to_vcf_line(
                    variants, variant_ids, max_alleles=max_alleles
                )
                if record is None:
                    variants = []
                    variants_ids = []
                    continue
                else:
                    self._append_to_clusters((record, variants, variant_ids))

                self._merge_last_records_in_clusters(max_alleles)
                self._print_first_clusters(f_main)

                variants = [var]
                variant_ids = [var_id]

            record = self._variant_cluster_to_vcf_line(
                variants, variant_ids, max_alleles=max_alleles
            )
            if record is not None:
                self._append_to_clusters((record, variants, variant_ids))
            self._merge_last_records_in_clusters(max_alleles)
            for c in self.clusters:
                print(c[0], file=f_main)
