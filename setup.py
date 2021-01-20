import glob
from setuptools import setup, find_packages


with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="cluster_vcf_records",
    version="0.13.1",
    description="Package to cluster VCF records. For use by gramtools and minos",
    packages=find_packages(),
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/cluster_vcf_records",
    test_suite="nose.collector",
    tests_require=["pytest"],
    install_requires=install_requires,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
