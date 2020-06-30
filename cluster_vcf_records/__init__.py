from pkg_resources import get_distribution

try:
    __version__ = get_distribution("cluster_vcf_records").version
except:
    __version__ = "local"


__all__ = [
    "allele_combinations",
    "utils",
    "variant_tracking",
    "vcf_clusterer",
    "vcf_file_read",
    "vcf_merge",
    "vcf_record",
    "vcf_record_cluster",
]

from cluster_vcf_records import *
from .vcf_clusterer import cluster
