import logging
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
    command = f"vcfbreakmulti {vcf_in} | vcfallelicprimitives -L 10000 | vt normalize -r {ref_fasta} - | vcfuniq > {vcf_out}"
    syscall(command)


def rm_rf(filename):
    syscall("rm -rf {filename}")
