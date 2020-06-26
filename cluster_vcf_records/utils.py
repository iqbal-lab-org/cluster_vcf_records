import logging
import os
import re
import subprocess
import sys

from cluster_vcf_records import vcf_record


def syscall(command):
    logging.debug("Run command: {command}")
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


def simplify_vcf(infile, outfile):
    """Removes records that are all null calls and ref calls.
    If record has GT, removes all non-called alleles and replaces FORMAT
    with just GT"""
    with open(infile) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                print(line, end="", file=f_out)
            else:
                record = vcf_record.VcfRecord(line)
                print("record:", record)
                if "GT" in record.FORMAT:
                    gt_indexes = set(re.split("[/|]", record.FORMAT["GT"]))
                    print("gt_indexes:", gt_indexes)
                    if "." in gt_indexes:
                        continue
                    gt_indexes = {int(x) for x in gt_indexes}
                    print("gt_indexes:", gt_indexes)
                    if gt_indexes == {0}:
                        continue

                    gt_indexes = sorted(list(gt_indexes))
                    record.ALT = [record.ALT[i - 1] for i in gt_indexes if i > 0]
                    record.FORMAT.clear()
                    if len(gt_indexes) == 1:
                        record.set_format_key_value("GT", "1/1")
                    else:
                        assert len(gt_indexes) == 2
                        if 0 in gt_indexes:
                            record.set_format_key_value("GT", "0/1")
                        else:
                            record.set_format_key_value("GT", "1/2")
                print(record, file=f_out)


def normalise_vcf(vcf_in, ref_fasta, vcf_out):
    # Would be nice to pipe all these to save disk IO. ie
    # f"vcfbreakmulti {vcf_in} | vcfallelicprimitives -L 10000 | vt normalize -r {ref_fasta} - | vcfuniq > {vcf_out}"
    # But am concerned about errors, where error code getting lost in pipes.
    # So run one by one using temp files and then clean up.
    vcf_breakmulti = f"{vcf_out}.1.breakmulti.vcf"
    vcf_allelic = f"{vcf_out}.2.allelicprimitives.vcf"
    vcf_normalize = f"{vcf_out}.3.normalize.vcf"
    syscall(f"vcfbreakmulti {vcf_in} > {vcf_breakmulti}")
    syscall(f"vcfallelicprimitives -L 10000 {vcf_breakmulti} > {vcf_allelic}")
    syscall(f"vt normalize -r {ref_fasta} {vcf_allelic} > {vcf_normalize}")
    syscall(f"vcfuniq {vcf_normalize} > {vcf_out}")
    os.unlink(vcf_breakmulti)
    os.unlink(vcf_allelic)
    os.unlink(vcf_normalize)


def rm_rf(filename):
    syscall(f"rm -rf {filename}")
