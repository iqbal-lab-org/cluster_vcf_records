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

from cluster_vcf_records import utils, vcf_file_read, vcf_record


Variant = namedtuple("Variant", ["seq_id", "pos", "ref", "alt"])


def _load_one_vcf_file(vcf_file, ref_seqs, ref_seq_to_id, ref_fasta, temp_dir):
    sample = vcf_file_read.get_sample_name_from_vcf_file(vcf_file)
    if sample is None:
        raise Exception(f"Error getting sample name from vcf file {vcf_file}")
    tmpdir = tempfile.mkdtemp(prefix="normalize_vcf.", dir=temp_dir)
    normalized_vcf = os.path.join(tmpdir, "normalized.vcf")
    utils.normalise_vcf(vcf_file, ref_fasta, normalized_vcf)
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
                ref_seq_index = ref_seq_to_id[record.CHROM]
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
        logging.info("Saving batch of data to disk {outfile}")
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


def load_slice_of_block(infile, seq_id, start, end):
    """Returns dictionary of variant id -> bitarray.
    start and end should be zero-based and include end points"""
    variants = {}
    if not isinstance(infile, pysam.libctabix.TabixFile):
        infile = pysam.TabixFile(infile)

    for line in infile.fetch(str(seq_id), start, end + 1):
        _, _, var_id, array = line.rstrip().split()
        variants[int(var_id)] = bitarray(array)
    return variants


class VariantTracker:
    def __init__(self, root_dir, ref_fasta, mem_limit=2):
        self.root_dir = os.path.abspath(root_dir)
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.mem_lmit = mem_limit
        self.variants = Variants()
        self.ref_seqs, self.ref_seq_to_id = VariantTracker.load_ref_seq_data(
            self.ref_fasta
        )
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
        return seqs, ref_seq_to_id

    def _load_data_from_root(self):
        with open(self.metadata_file) as f:
            self.var_block_files, self.samples = json.load(f)
        self.variants.load_from_file(self.variants_file)

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

    def construct_from_vcf_files(
        self, infiles, temp_dir, cpus=1, mem_limit=2, force=False, sample_limit=None
    ):
        if force:
            utils.rm_rf(self.root_dir)
        os.mkdir(self.root_dir)
        self.var_block_files = []
        self.samples = []
        self._add_var_block_file()
        self.variants = Variants()
        var_block = VariantBlock()

        for i in range(0, len(infiles), cpus):
            with multiprocessing.Pool(cpus) as pool:
                new_variants_lists = pool.starmap(
                    _load_one_vcf_file,
                    zip(
                        infiles[i : i + cpus],
                        itertools.repeat(self.ref_seqs),
                        itertools.repeat(self.ref_seq_to_id),
                        itertools.repeat(self.ref_fasta),
                        itertools.repeat(temp_dir),
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

        logging.info(f"Loaded all {len(infiles)} VCF files")
        self._write_last_var_block_file(var_block)

        logging.info(f"Saving metadata to file {self.metadata_file}")
        self._save_metadata_file()

        logging.info(f"Saving variants to file {self.variants_file}")
        self.variants.save_to_file(self.variants_file)
        logging.info("Finished")

    def cluster(self, outdir):
        # FIXME
        pass
