import logging
import os
import subprocess
import sys

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
