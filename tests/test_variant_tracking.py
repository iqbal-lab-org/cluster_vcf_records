import filecmp
import os
import pytest

from bitarray import bitarray
import pyfastaq

from cluster_vcf_records import utils, variant_tracking, vcf_file_read, vcf_record

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "variant_tracking")


def test_vcf_records_make_same_allele_combination():
    ref_seqs = {"ref1": "GCTGT"}
    record1 = vcf_record.VcfRecord("ref1\t1\t.\tGCT\tGC,GCGT\t.\t.\t.")
    record2 = vcf_record.VcfRecord("ref1\t5\t.\tT\tTGG,G\t.\t.\t.")
    record3 = vcf_record.VcfRecord("ref2\t5\t.\tT\tTGG,G\t.\t.\t.")
    assert variant_tracking.vcf_records_make_same_allele_combination(
        record1, record2, ref_seqs
    )
    assert not variant_tracking.vcf_records_make_same_allele_combination(
        record1, record3, ref_seqs
    )


def test_variants_overlap():
    v1 = variant_tracking.Variant(0, 10, "A", "T")
    v2 = variant_tracking.Variant(0, 10, "AG", "T")
    v3 = variant_tracking.Variant(1, 10, "A", "T")
    v4 = variant_tracking.Variant(0, 11, "A", "T")
    v5 = variant_tracking.Variant(0, 12, "A", "T")
    f = variant_tracking.variants_overlap
    assert f(v1, v2)
    assert f(v2, v1)
    assert not f(v1, v3)
    assert not f(v3, v1)
    assert not f(v1, v4)
    assert not f(v4, v1)
    assert not f(v1, v5)
    assert not f(v5, v1)
    assert f(v2, v4)
    assert f(v4, v2)
    assert not f(v2, v5)
    assert not f(v5, v2)


def test_Variants():
    variants = variant_tracking.Variants()
    assert len(variants) == 0
    variant1 = variant_tracking.Variant(1, 3, "A", "T")
    assert variant1 not in variants
    var_id = variants.add(variant1)
    assert variant1 in variants
    assert variants[variant1] == var_id
    variant2 = variant_tracking.Variant(1, 2, "A", "T")
    variant3 = variant_tracking.Variant(0, 1, "A", "T")
    variants.add(variant2)
    variants.add(variant3)
    expect = [
        (variant3, variants[variant3]),
        (variant2, variants[variant2]),
        (variant1, variants[variant1]),
    ]
    assert list(variants.sorted_iter()) == expect

    tmp_file = "tmp.variants.save_to_file"
    utils.rm_rf(tmp_file)
    variants.save_to_file(tmp_file)
    new_variants = variant_tracking.Variants()
    new_variants.load_from_file(tmp_file)
    assert variants == new_variants
    os.unlink(tmp_file)


def test_VariantBlock():
    # Contruct and add a variant and sample
    block = variant_tracking.VariantBlock()
    assert block.number_of_samples() == 0
    assert block.number_of_variants() == 0
    block.add_variants(1)
    assert block.number_of_samples() == 0
    assert block.number_of_variants() == 1
    block.add_samples(1)
    assert block.number_of_samples() == 1
    assert block.number_of_variants() == 1

    # Getting and setting variant
    assert not block.has_variant(0, 0)
    block.set_variant(0, 0)
    assert block.has_variant(0, 0)

    # Add more samples and variant, check first variant and sample not changed
    block.add_variants(3)
    block.add_samples(2)
    assert block.number_of_samples() == 3
    assert block.number_of_variants() == 4
    assert block.has_variant(0, 0)
    assert not block.has_variant(0, 1)
    assert not block.has_variant(0, 2)
    assert not block.has_variant(1, 0)
    assert not block.has_variant(1, 1)
    assert not block.has_variant(1, 2)
    assert not block.has_variant(2, 0)
    assert not block.has_variant(2, 1)
    assert not block.has_variant(2, 2)
    assert not block.has_variant(3, 0)
    assert not block.has_variant(3, 1)
    assert not block.has_variant(3, 2)
    block.set_variant(2, 2)
    block.set_variant(3, 0)
    block.set_variant(3, 1)

    # Save to file
    variants = variant_tracking.Variants()
    variants.add(variant_tracking.Variant(0, 0, "A", "G"))
    variants.add(variant_tracking.Variant(0, 2, "G", "T"))
    variants.add(variant_tracking.Variant(0, 2, "G", "C"))
    variants.add(variant_tracking.Variant(1, 42, "G", "C"))
    tmp_file = "tmp.variant_tracking.block.tsv.gz"
    utils.rm_rf(tmp_file)
    utils.rm_rf(tmp_file + ".tbi")
    block.write_to_bgzip_file_and_tab_index(tmp_file, variants)
    wanted_ids = set([v for k, v in variants.sorted_iter()])

    # Load slices from file. Note that none of the variants had variants[1], so
    # should not be in the file
    assert variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 1, 0, 0) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 1, 41, 41) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 1, 43, 43) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 0, 1, 1) == {}
    expect_vars = {0: block.bitarrays[0]}
    assert (
        variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 0, 0, 0)
        == expect_vars
    )
    assert variant_tracking.load_slice_of_block(tmp_file, {1}, 0, 0, 0) == {}
    assert (
        variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 0, 0, 1)
        == expect_vars
    )
    expect_vars[2] = bitarray(block.bitarrays[2])
    assert (
        variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 0, 0, 2)
        == expect_vars
    )
    assert (
        variant_tracking.load_slice_of_block(tmp_file, wanted_ids, 0, 0, 3)
        == expect_vars
    )

    # Load variant patterns from slices of block. Make another block file to test
    # getting from >1 file
    block.clear_samples()
    variants.add(variant_tracking.Variant(0, 1, "C", "G"))
    variants.add(variant_tracking.Variant(0, 10, "T", "A"))
    block.add_variants(2)
    block.add_samples(2)
    block.set_variant(0, 1)
    block.set_variant(1, 0)
    block.set_variant(1, 1)
    block.set_variant(4, 0)
    block.set_variant(5, 1)
    tmp_file2 = "tmp.variant_tracking.block.2.tsv.gz"
    utils.rm_rf(tmp_file2)
    utils.rm_rf(tmp_file2 + ".tbi")
    block.write_to_bgzip_file_and_tab_index(tmp_file2, variants)
    wanted_ids = set([v for k, v in variants.sorted_iter()])
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 1, 0, 41
    )
    assert got_patterns == set()
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 1, 0, 42
    )
    assert got_patterns == {(3,)}
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 1, 42, 42
    )
    assert got_patterns == {(3,)}
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 1, 42, 43
    )
    assert got_patterns == {(3,)}
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 1, 43, 43
    )
    assert got_patterns == set()
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 0, 0, 9
    )
    expect_patterns = {(0, 1), (2,), (0,), (1, 4)}
    assert got_patterns == expect_patterns
    got_patterns = variant_tracking.var_patterns_from_block_slices(
        [tmp_file, tmp_file2], wanted_ids, 0, 0, 10
    )
    expect_patterns = {(0, 1, 5), (2,), (0,), (1, 4)}
    assert got_patterns == expect_patterns

    os.unlink(tmp_file)
    os.unlink(tmp_file + ".tbi")
    os.unlink(tmp_file2)
    os.unlink(tmp_file2 + ".tbi")


def test_var_pattern_to_alleles():
    ref_seq = pyfastaq.sequences.Fasta("ref", "GTAGCTGACTGCTAGTGTA")
    #                                          0123456789012345678
    variants = {
        0: variant_tracking.Variant(0, 3, "G", "T"),
        1: variant_tracking.Variant(0, 3, "G", "A"),
        2: variant_tracking.Variant(0, 4, "C", "A"),
        3: variant_tracking.Variant(0, 6, "GA", "G"),
        4: variant_tracking.Variant(0, 8, "C", "CAAA"),
        5: variant_tracking.Variant(0, 9, "T", "CG"),
        6: variant_tracking.Variant(0, 10, "G", "C"),
        7: variant_tracking.Variant(0, 10, "GC", "G"),
    }
    f = variant_tracking.var_pattern_to_alleles
    assert f(variants, {0}, ref_seq, 2, 3) == {"AT"}
    assert f(variants, {0}, ref_seq, 2, 4) == {"ATC"}
    assert f(variants, {0}, ref_seq, 3, 3) == {"T"}
    assert f(variants, {0}, ref_seq, 3, 4) == {"TC"}
    assert f(variants, {0, 1}, ref_seq, 3, 4) == {"TC", "AC"}
    assert f(variants, {0, 2}, ref_seq, 3, 4) == {"TA"}
    assert f(variants, {0, 2, 3}, ref_seq, 3, 7) == {"TATG"}
    assert f(variants, {0, 2, 3}, ref_seq, 3, 8) == {"TATGC"}
    assert f(variants, {0, 2, 3}, ref_seq, 3, 9) == {"TATGCT"}
    assert f(variants, {0, 2, 3, 4}, ref_seq, 3, 9) == {"TATGCAAAT"}
    assert f(variants, {0, 2, 3, 4, 5}, ref_seq, 3, 9) == {"TATGCAAACG"}
    assert f(variants, {0, 2, 3, 4, 5}, ref_seq, 3, 10) == {"TATGCAAACGG"}
    assert f(variants, {0, 2, 3, 4, 5, 6}, ref_seq, 3, 10) == {"TATGCAAACGC"}
    assert f(variants, {0, 2, 3, 4, 5, 6}, ref_seq, 3, 11) == {"TATGCAAACGCC"}
    assert f(variants, {0, 2, 3, 4, 5, 7}, ref_seq, 3, 11) == {"TATGCAAACGG"}


def test_load_one_vcf_file():
    vcf_file = os.path.join(data_dir, "load_one_vcf_file.vcf")
    ref_fasta = os.path.join(data_dir, "load_one_vcf_file.fa")
    (
        ref_seqs,
        ref_names,
        ref_seq_to_id,
    ) = variant_tracking.VariantTracker.load_ref_seq_data(ref_fasta)
    tmp_dir = "tmp.load_one_vcf_file"
    utils.rm_rf(tmp_dir)
    os.mkdir(tmp_dir)
    got_sample, got_variants = variant_tracking._load_one_vcf_file(
        vcf_file, ref_seqs, ref_seq_to_id, ref_fasta, tmp_dir, True
    )
    assert got_sample == "sample_42"
    expect_variants = [
        variant_tracking.Variant(seq_id=0, pos=1, ref="T", alt="TCGC"),
        variant_tracking.Variant(seq_id=0, pos=2, ref="C", alt="G"),
        variant_tracking.Variant(seq_id=0, pos=6, ref="A", alt="T"),
        variant_tracking.Variant(seq_id=0, pos=8, ref="T", alt="G"),
        variant_tracking.Variant(seq_id=1, pos=1, ref="G", alt="C"),
        variant_tracking.Variant(seq_id=1, pos=1, ref="G", alt="A"),
    ]
    assert got_variants == expect_variants
    os.rmdir(tmp_dir)


def test_VariantTracker_make_from_vcf_then_save_then_load_then_cluster():
    # Create from VCF files and save to directory
    vcf_files = [
        os.path.join(data_dir, f"construct_from_vcf_files.{i}.vcf") for i in (1, 2, 3)
    ]
    ref_fasta = os.path.join(data_dir, "construct_from_vcf_files.ref.fa")
    root_dir = "tmp.construct_from_vcf_files"
    utils.rm_rf(root_dir)
    tmp_dir = f"{root_dir}.tmp"
    utils.rm_rf(tmp_dir)
    os.mkdir(tmp_dir)
    root_dir = "tmp.construct_from_vcf_files"
    tracker = variant_tracking.VariantTracker(root_dir, ref_fasta)
    tracker.merge_vcf_files(vcf_files, temp_dir=tmp_dir, cpus=2, sample_limit=2)

    # Check variables ok and files made ok
    expect_samples = [["sample_1", "sample_2"], ["sample_3"]]
    expect_block_files = ["block.0.tsv.gz", "block.1.tsv.gz"]
    assert tracker.samples == expect_samples
    assert tracker.var_block_files == expect_block_files
    expect_dir = os.path.join(data_dir, "construct_from_vcf_files.expect")
    expect_variants = variant_tracking.Variants()
    expect_variants.load_from_file(os.path.join(expect_dir, "variants.tsv.gz"))
    assert tracker.variants == expect_variants
    for block_file in expect_block_files:
        expect = os.path.join(expect_dir, block_file)
        got = os.path.join(root_dir, block_file)
        with vcf_file_read.open_vcf_file_for_reading(expect) as f:
            expect_lines = [x.rstrip() for x in f]
        with vcf_file_read.open_vcf_file_for_reading(got) as f:
            got_lines = [x.rstrip() for x in f]
        assert expect_lines == got_lines
        assert os.path.exists(f"{got}.tbi")

    # Now the root_dir exists, constructor should load data from directory
    tracker = variant_tracking.VariantTracker(root_dir, ref_fasta)
    assert tracker.samples == expect_samples
    assert tracker.var_block_files == expect_block_files
    assert tracker.variants == expect_variants

    # Run the clustering. max_ref_len 8 is length of longest REF string in
    # the test data
    cluster_out = "tmp.variant_tracker.cluster"
    vcf_out = f"{cluster_out}.vcf"
    excluded_out = f"{cluster_out}.excluded.tsv"
    tracker.cluster(cluster_out, 10)
    expect_vcf = os.path.join(data_dir, "cluster.max_ref_8.vcf")
    expect_excluded = os.path.join(data_dir, "cluster.max_ref_8.excluded.tsv")
    assert filecmp.cmp(vcf_out, expect_vcf, shallow=False)
    assert filecmp.cmp(excluded_out, expect_excluded, shallow=False)
    os.unlink(vcf_out)
    os.unlink(excluded_out)

    # Set ref length limit to exclude the ref length 8 variant
    tracker.cluster(cluster_out, 7)
    expect_vcf = os.path.join(data_dir, "cluster.max_ref_7.vcf")
    expect_excluded = os.path.join(data_dir, "cluster.max_ref_7.excluded.tsv")
    assert filecmp.cmp(vcf_out, expect_vcf, shallow=False)
    assert filecmp.cmp(excluded_out, expect_excluded, shallow=False)
    os.unlink(vcf_out)
    os.unlink(excluded_out)

    # Set allele limit to trigger on the length 8 variant, so we only
    # get the alt alleles for each sample
    tracker.cluster(cluster_out, 8, max_alleles=5)
    expect_vcf = os.path.join(data_dir, "cluster.max_alleles_5.vcf")
    expect_excluded = os.path.join(data_dir, "cluster.max_alleles_5.excluded.tsv")
    assert filecmp.cmp(vcf_out, expect_vcf, shallow=False)
    assert filecmp.cmp(excluded_out, expect_excluded, shallow=False)
    os.unlink(vcf_out)
    os.unlink(excluded_out)

    utils.rm_rf(tmp_dir)
    utils.rm_rf(root_dir)
