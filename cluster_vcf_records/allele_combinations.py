import copy
import itertools
from operator import attrgetter

from cluster_vcf_records import variant_tracking


def var_cluster_to_coords_and_snps_and_non_snps(variants):
    assert len(variants) > 0
    nucleotides = {"A", "C", "G", "T"}
    snps = {}  # position => set of alts (and the ref nucleotide)
    non_snps = []  #  list of variants
    start = float("Inf")
    end = -1

    for var in variants:
        start = min(start, var.pos)
        end = max(end, var.pos + len(var.ref) - 1)

        if var.ref in nucleotides and var.alt in nucleotides:
            if var.pos not in snps:
                snps[var.pos] = {var.ref}
            snps[var.pos].add(var.alt)
        else:
            ref = var.ref
            alt = var.alt
            i = 0
            while i < len(ref) - 1 and i < len(alt) and ref[i] == alt[i]:
                i += 1
            non_snps.append(
                variant_tracking.Variant(var.seq_id, var.pos + i, ref[i:], alt[i:])
            )

    return start, end, snps, non_snps


def any_vars_overlap(variants):
    if len(variants) < 2:
        return False

    for v1, v2 in itertools.combinations(variants, 2):
        if variant_tracking.variants_overlap(v1, v2):
            return True

    return False


def var_cluster_to_coords_and_alts(variants, ref_seq, max_alleles=None):
    if len(variants) == 1:
        end = variants[0].pos + len(variants[0].ref) - 1
        return variants[0].pos, end, {variants[0].alt}

    if max_alleles is None:
        max_alleles = float("Inf")
    (
        final_start,
        final_end,
        snps,
        non_snps,
    ) = var_cluster_to_coords_and_snps_and_non_snps(variants)

    # generate all the allele combinations from the SNPs.
    snp_positions = []
    snp_nucleotides = []
    for position in sorted(snps):
        snp_positions.append(position)
        snp_nucleotides.append(sorted(list(snps[position])))

    # work out min total alleles without making them. Making them could
    # take a long time if too many! Can only put lower bound on
    # the final unique number because all the combinations may have duplicates.
    total_alleles_lower_bound = 1
    for x in snp_nucleotides:
        total_alleles_lower_bound *= len(x)
    if len(non_snps):
        total_alleles_lower_bound *= len(non_snps)

    if total_alleles_lower_bound > max_alleles:
        return final_start, final_end, None

    alts = set()
    ref_seq_for_vcf = ref_seq[final_start : final_end + 1]

    # This generates all possible combinations of indels. Checks if
    # each combination contains any that overlap, in which case that combination
    # is not used because can't apply two conflicting indels.
    overlaps = {}
    for v1, v2 in itertools.combinations(non_snps, 2):
        if variant_tracking.variants_overlap(v1, v2):
            key1, key2 = sorted([v1, v2])
            if key1 not in overlaps:
                overlaps[key1] = set()
            overlaps[key1].add(key2)

    # snp_nucleotides has the ref and all the alt snps.
    # This loops over all combinations of them, so we're getting every
    # possibility of using each snp or not.
    for combination in itertools.product(*snp_nucleotides):
        alt_seq = list(ref_seq_for_vcf)
        for i, position in enumerate(snp_positions):
            alt_seq[position - final_start] = combination[i]
        alts.add("".join(alt_seq))
        if len(alts) > max_alleles:
            return final_start, final_end, None

        number_of_combos_used = 1
        for number_of_non_snps in range(1, len(non_snps) + 1):
            # If the previous combination of number_of_non_snps of non_snps
            # resulted in none getting used because of overlapping non_snps,
            # then taking more than number_of_non_snps won't work either so stop
            if number_of_combos_used == 0:
                break
            else:
                number_of_combos_used = 0

            for non_snp_combo in itertools.combinations(non_snps, number_of_non_snps):
                overlap = False
                if number_of_non_snps > 1:
                    for v1 in non_snps:
                        if v1 in overlaps:
                            for v2 in non_snps:
                                if v2 in overlaps[v1]:
                                    overlap = True
                                    break

                        if overlap:
                            break

                if overlap:
                    continue

                number_of_combos_used += 1
                non_snp_combo = sorted(list(non_snp_combo), key=attrgetter("pos"))
                new_alt = copy.copy(alt_seq)
                for non_snp in reversed(non_snp_combo):
                    start_pos = non_snp.pos - final_start
                    new_alt[start_pos : start_pos + len(non_snp.ref)] = non_snp.alt
                    alts.add("".join(new_alt))
                    if len(alts) > max_alleles:
                        return final_start, final_end, None

    alts.discard(ref_seq_for_vcf)
    return final_start, final_end, alts
