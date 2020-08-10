import logging
import os
import re
import shutil
import subprocess
import sys

import pyfastaq

from cluster_vcf_records import vcf_record, vcf_file_read


def syscall(command):
    logging.debug(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "\nOutput from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "\nOutput from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    return completed_process


def simplify_vcf(
    infile, outfile, ref_seqs=None, keep_ref_calls=False, max_ref_call_alleles=10
):
    """Removes records that are all null calls and (optionally) ref calls.
    If record has GT, removes all non-called alleles and replaces FORMAT
    with just GT.
    If ref_Seqs is given, should be a dictionary of reference name -> sequence.
    Any calls where the REF string does not match the sequence is removed."""
    f_in = vcf_file_read.open_vcf_file_for_reading(infile)
    with open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                print(line, end="", file=f_out)
            else:
                try:
                    record = vcf_record.VcfRecord(line)
                except:
                    logging.warning(
                        f"Bad VCF line from file {infile}. Ignoring: {line.rstrip()}"
                    )
                    continue

                if ref_seqs is not None and not (
                    record.ref_string_matches_dict_of_ref_sequences(ref_seqs)
                ):
                    continue

                if "GT" in record.FORMAT:
                    gt_indexes = set(re.split("[/|]", record.FORMAT["GT"]))
                    if "." in gt_indexes:
                        continue
                    gt_indexes = {int(x) for x in gt_indexes}
                    record.FORMAT.clear()
                    if gt_indexes == {0}:
                        if keep_ref_calls and len(record.ALT) <= max_ref_call_alleles:
                            record.set_format_key_value("GT", "0/0")
                        else:
                            continue
                    else:
                        gt_indexes = sorted(list(gt_indexes))
                        record.ALT = [record.ALT[i - 1] for i in gt_indexes if i > 0]
                        if len(gt_indexes) == 1:
                            record.set_format_key_value("GT", "1/1")
                        else:
                            assert len(gt_indexes) == 2
                            if 0 in gt_indexes:
                                record.set_format_key_value("GT", "0/1")
                            else:
                                record.set_format_key_value("GT", "1/2")
                print(record, file=f_out)

    f_in.close()


def normalise_vcf(vcf_in, ref_fasta, vcf_out, break_alleles=True):
    # Would be nice to pipe all these to save disk IO. ie
    # f"vcfbreakmulti {vcf_in} | vcfallelicprimitives -L 10000 | vt normalize -r {ref_fasta} - | vcfuniq > {vcf_out}"
    # But am concerned about errors, where error code getting lost in pipes.
    # So run one by one using temp files and then clean up.
    vcf_breakmulti = f"{vcf_out}.1.breakmulti.vcf"
    vcf_allelic = f"{vcf_out}.2.allelicprimitives.vcf"
    vcf_normalize = f"{vcf_out}.3.normalize.vcf"
    syscall(f"vcfbreakmulti {vcf_in} > {vcf_breakmulti}")
    if break_alleles:
        syscall(f"vcfallelicprimitives -L 10000 {vcf_breakmulti} > {vcf_allelic}")
    else:
        vcf_allelic = vcf_breakmulti
    syscall(f"vt normalize -r {ref_fasta} {vcf_allelic} > {vcf_normalize}")
    syscall(f"vcfuniq {vcf_normalize} > {vcf_out}")
    os.unlink(vcf_breakmulti)
    if break_alleles:
        os.unlink(vcf_allelic)
    os.unlink(vcf_normalize)


def rm_rf(filename):
    syscall(f"rm -rf {filename}")


def cat_vcfs(vcf_files, outfile, delete_files=False):
    shutil.copyfile(vcf_files[0], outfile)
    if delete_files:
        os.unlink(vcf_files[0])
    with open(outfile, "a") as f_out:
        for vcf in vcf_files[1:]:
            with open(vcf) as f_in:
                for line in  f_in:
                    if not line.startswith("#"):
                        print(line, end="", file=f_out)
            if delete_files:
                os.unlink(vcf)


def cat_tsvs(tsv_files, outfile, delete_files=False):
    shutil.copyfile(tsv_files[0], outfile)
    if delete_files:
        os.unlink(tsv_files[0])
    with open(outfile, "a") as f_out:
        for tsv in tsv_files[1:]:
            with open(tsv) as f_in:
                for i, line in  enumerate(f_in):
                    if i > 1:
                        print(line, end="", file=f_out)
            if delete_files:
                os.unlink(tsv)
