import glob
from setuptools import setup, find_packages


setup(
    name="cluster_vcf_records",
    version="0.10.2",
    description="Package to cluster VCF records. For use by gramtools and minos",
    packages=find_packages(),
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/cluster_vcf_records",
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=["pyfastaq >= 3.14.0"],
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
