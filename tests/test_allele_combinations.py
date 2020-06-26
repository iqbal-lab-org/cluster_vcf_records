import pytest

import pyfastaq

from cluster_vcf_records import allele_combinations, variant_tracking

def test_var_cluster_to_coords_and_snps_and_non_snps():
    variants = [variant_tracking.Variant(0, 10, "A", "G")]
    f = allele_combinations.var_cluster_to_coords_and_snps_and_non_snps
    snps = {10: {"A", "G"}}
    assert f(variants) == (10, 10, snps, [])

    indel = variant_tracking.Variant(0, 10, "AT", "A")
    alt_indel = variant_tracking.Variant(0, 11, "T", "")
    variants.extend([
        variant_tracking.Variant(0, 10, "A", "T"),
        indel,
    ])
    snps[10].add("T")
    assert f(variants) == (10, 11, snps, [alt_indel])

def test_any_vars_overlap():
    variants = [variant_tracking.Variant(0, 4, "T", "G")]
    assert not allele_combinations.any_vars_overlap(variants)
    variants.append(variant_tracking.Variant(0, 5, "TA", "G"))
    assert not allele_combinations.any_vars_overlap(variants)
    variants.append(variant_tracking.Variant(0, 6, "A", "G"))
    assert  allele_combinations.any_vars_overlap(variants)


def test_var_cluster_to_coords_and_alts():
    ref_seq = pyfastaq.sequences.Fasta("ref", "CTAGTCGATGCACTGATAGTA")
    #                                          012345678901234567890
    variants = [variant_tracking.Variant(0, 4, "T", "G")]
    f = allele_combinations.var_cluster_to_coords_and_alts
    assert f(variants, ref_seq) == (4, 4, {"G"})

    variants.append(variant_tracking.Variant(0, 4, "T", "A"))
    assert f(variants, ref_seq) == (4, 4, {"A", "G"})

    variants = [
        variant_tracking.Variant(0, 4, "T", "G"),
        variant_tracking.Variant(0, 5, "C", "A"),
    ]
    assert f(variants, ref_seq) == (4, 5, {"TA", "GC", "GA"})

    variants = [
        variant_tracking.Variant(0, 4, "TC", "T"),
        variant_tracking.Variant(0, 5, "C", "A"),
    ]
    assert f(variants, ref_seq) == (4, 5, {"TA", "T"})

    variants = [
        variant_tracking.Variant(0, 4, "TC", "T"),
        variant_tracking.Variant(0, 4, "TCG", "T"),
        variant_tracking.Variant(0, 5, "C", "A"),
    ]
    assert f(variants, ref_seq) == (4, 6, {"TG", "T", "TAG"})

    # Annoying edge case. The two deletions appear to overlap, because
    # we're using the VCF convention of using the nucleotide before the deletion
    # in the ref and alt. Check that we get the result of both deletions
    # being applied
    variants = [
        variant_tracking.Variant(0, 4, "TCG", "T"),
        variant_tracking.Variant(0, 6, "GAT", "G"),
    ]
    assert f(variants, ref_seq) == (4, 8, {"TAT", "TCG", "T"})

    variants.append(variant_tracking.Variant(0, 5, "C", "A"))
    assert f(variants, ref_seq) == (4, 8, {"TAT", "TCG", "T", "TAG", "TAGAT"})
