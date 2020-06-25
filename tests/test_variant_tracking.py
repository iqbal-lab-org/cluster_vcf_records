import filecmp
import os
import pytest

from bitarray import bitarray
import pyfastaq

from cluster_vcf_records import utils, variant_tracking

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "variant_tracking")


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
    tmp_file = "tmp.variant_trcking.block.tsv.gz"
    utils.rm_rf(tmp_file)
    utils.rm_rf(tmp_file + ".tbi")
    block.write_to_bgzip_file_and_tab_index(tmp_file, variants)

    # Load slices from file. Note that none of the variants had variants[1], so
    # should not be in the file
    assert variant_tracking.load_slice_of_block(tmp_file, 1, 0, 0) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, 1, 41, 41) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, 1, 43, 43) == {}
    assert variant_tracking.load_slice_of_block(tmp_file, 0, 1, 1) == {}
    expect_vars = {0: block.bitarrays[0]}
    assert variant_tracking.load_slice_of_block(tmp_file, 0, 0, 0) == expect_vars
    assert variant_tracking.load_slice_of_block(tmp_file, 0, 0, 1) == expect_vars
    expect_vars[2] = bitarray(block.bitarrays[2])
    assert variant_tracking.load_slice_of_block(tmp_file, 0, 0, 2) == expect_vars
    assert variant_tracking.load_slice_of_block(tmp_file, 0, 0, 3) == expect_vars

    os.unlink(tmp_file)
    os.unlink(tmp_file + ".tbi")


def test_load_one_vcf_file():
    vcf_file = os.path.join(data_dir, "load_one_vcf_file.vcf")
    ref_fasta = os.path.join(data_dir, "load_one_vcf_file.fa")
    ref_seqs, ref_seq_to_id = variant_tracking.VariantTracker.load_ref_seq_data(
        ref_fasta
    )
    tmp_dir = "tmp.load_one_vcf_file"
    utils.rm_rf(tmp_dir)
    os.mkdir(tmp_dir)
    got_sample, got_variants = variant_tracking._load_one_vcf_file(
        vcf_file, ref_seqs, ref_seq_to_id, ref_fasta, tmp_dir
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


def test_VariantTracker_make_from_vcf_then_save_then_load():
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
    tracker.construct_from_vcf_files(vcf_files, tmp_dir, cpus=2, sample_limit=2)

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
        assert filecmp.cmp(got, expect, shallow=False)
        assert os.path.exists(f"{got}.tbi")

    # Now the root_dir exists, constructor should load data from directory
    tracker = variant_tracking.VariantTracker(root_dir, ref_fasta)
    assert tracker.samples == expect_samples
    assert tracker.var_block_files == expect_block_files
    assert tracker.variants == expect_variants
    utils.rm_rf(tmp_dir)
    utils.rm_rf(root_dir)
