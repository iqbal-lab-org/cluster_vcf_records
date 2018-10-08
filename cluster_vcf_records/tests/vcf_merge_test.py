import datetime
import filecmp
import os
import unittest

from cluster_vcf_records import vcf_file_read, vcf_merge
from cluster_vcf_records import __version__ as cluster_vcf_records_version

modules_dir = os.path.dirname(os.path.abspath(vcf_merge.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'vcf_merge')


def check_vcfs(expected_vcf, got_vcf):
    expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(expected_vcf)
    got_header, got_vcf_records = vcf_file_read.vcf_file_to_list(got_vcf)
    for i in range(len(expected_header)):
        if expected_header[i].startswith('##fileDate='):
            expected_header[i] = '##fileDate=' + str(datetime.date.today())
        elif expected_header[i].startswith('##source=minos'):
            expected_header[i] = '##source=cluster_vcf_records, version ' + cluster_vcf_records_version

    return expected_header == got_header and expected_vcf_records == got_vcf_records


class TestVcfMerger(unittest.TestCase):
    def test_dict_of_vars_to_vcf_file(self):
        '''test _dict_of_vars_to_vcf_file'''
        variants = {
            'ref_42': {
                10: {'A': {'G', 'T', 'TT'}, 'AC': {'G'}},
                11: {'C': {'G'}},
            },
            'ref_43': {
                41: {'T': {'G'}},
            },
        }
        tmp_file = 'tmp.test.dict_of_vars_to_vcf_file.out.vcf'
        vcf_merge._dict_of_vars_to_vcf_file(variants, tmp_file)
        expected_file = os.path.join(data_dir, 'dict_of_vars_to_vcf_file.vcf')
        self.assertTrue(check_vcfs(expected_file, tmp_file))
        os.unlink(tmp_file)

