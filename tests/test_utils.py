import filecmp
import os
import pytest

import pyfastaq
from cluster_vcf_records import vcf_file_read, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_simplify_vcf():
    infile = os.path.join(data_dir, "simplify_vcf.in.vcf")
    ref_fa = os.path.join(data_dir, "simplify_vcf.ref.fa")
    ref_seqs = {}
    pyfastaq.tasks.file_to_dict(ref_fa, ref_seqs)
    tmp_out = "tmp.simplify_vcf.out.vcf"
    utils.rm_rf(tmp_out)
    utils.simplify_vcf(infile, tmp_out, ref_seqs=ref_seqs)
    expect = os.path.join(data_dir, "simplify_vcf.expect.vcf")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)
    utils.simplify_vcf(infile + ".gz", tmp_out, ref_seqs=ref_seqs)
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    utils.simplify_vcf(infile, tmp_out, keep_ref_calls=True, ref_seqs=ref_seqs)
    expect = os.path.join(data_dir, "simplify_vcf.expect_keep_ref_calls.vcf")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_normalise_vcf():
    infile = os.path.join(data_dir, "normalise_vcf.in.vcf")
    ref_fa = os.path.join(data_dir, "normalise_vcf.in.fa")
    expect = os.path.join(data_dir, "normalise_vcf.out.vcf")
    tmp_out = "tmp.normalise_vcf.vcf"
    utils.rm_rf(tmp_out)
    utils.normalise_vcf(infile, ref_fa, tmp_out)
    expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(expect)
    got_header, got_vcf_records = vcf_file_read.vcf_file_to_list(tmp_out)
    # The normalizing commands add lots of lines to the header.
    # We don't care about those, so just check the actual records.
    assert got_vcf_records == expected_vcf_records
    os.unlink(tmp_out)
