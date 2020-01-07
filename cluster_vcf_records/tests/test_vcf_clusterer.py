import filecmp
import os
import unittest

from cluster_vcf_records import vcf_clusterer, vcf_record, vcf_record_cluster

modules_dir = os.path.dirname(os.path.abspath(vcf_clusterer.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "vcf_clusterer")


class TestVcfClusterer(unittest.TestCase):
    def test_load_vcf_files(self):
        """test _load_vcf_files"""
        vcf_file_1 = os.path.join(data_dir, "load_vcf_files.1.vcf")

        expected_headers = {vcf_file_1: ["#file1 header1", "#file1 header2"]}

        expected_records = {
            "ref.1": [
                vcf_record.VcfRecord("ref.1\t5\tid3\tG\tA\tPASS\tSVTYPE=SNP\tGT\t1/1"),
                vcf_record.VcfRecord("ref.1\t10\tid1\tA\tT\tPASS\tSVTYPE=SNP\tGT\t1/1"),
            ],
            "ref.2": [
                vcf_record.VcfRecord("ref.2\t42\tid2\tG\tC\tPASS\tSVTYPE=SNP\tGT\t1/1")
            ],
        }

        expected_sample = "sample"
        got_sample, got_headers, got_records = vcf_clusterer.VcfClusterer._load_vcf_files(
            [vcf_file_1], None
        )
        self.assertEqual(expected_sample, got_sample)
        self.assertEqual(expected_headers, got_headers)
        self.assertEqual(expected_records, got_records)

        vcf_file_2 = os.path.join(data_dir, "load_vcf_files.2.vcf")
        expected_headers[vcf_file_2] = [
            "#file2 header",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_from_vcf_2",
        ]
        expected_records["ref.3"] = [
            vcf_record.VcfRecord("ref.3\t8\tid5\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1")
        ]
        expected_records["ref.1"].insert(
            1, vcf_record.VcfRecord("ref.1\t8\tid4\tC\tG\tPASS\tSVTYPE=SNP\tGT\t1/1")
        )
        expected_sample = "sample_from_vcf_2"
        got_sample, got_headers, got_records = vcf_clusterer.VcfClusterer._load_vcf_files(
            [vcf_file_1, vcf_file_2], None
        )
        self.assertEqual(expected_sample, got_sample)
        self.assertEqual(expected_headers, got_headers)
        self.assertEqual(expected_records, got_records)

    def test_expand_alts_in_vcf_record_list(self):
        """test _expand_alts_in_vcf_record_list"""
        vcf_records = [
            vcf_record.VcfRecord("ref\t42\t.\tA\tC,T\t.\tPASS\tSVTYPE=SNP\n"),
            vcf_record.VcfRecord("ref\t42\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n"),
            vcf_record.VcfRecord("ref\t44\t.\tA\tC\t.\tPASS\tSVTYPE=SNP\n"),
        ]
        expected = [
            vcf_record.VcfRecord("ref\t42\t.\tA\tC\t.\tPASS\tSVTYPE=SNP\n"),
            vcf_record.VcfRecord("ref\t42\t.\tA\tT\t.\tPASS\tSVTYPE=SNP\n"),
            vcf_record.VcfRecord("ref\t42\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n"),
            vcf_record.VcfRecord("ref\t44\t.\tA\tC\t.\tPASS\tSVTYPE=SNP\n"),
        ]

        got = vcf_clusterer.VcfClusterer._expand_alts_in_vcf_record_list(vcf_records)
        self.assertEqual(expected, got)

    def test_expand_alts_and_remove_duplicates_in_list(self):
        """test _expand_alts_and_remove_duplicates_in_list"""
        #          123456789012345678901234567890
        ref_seq = "ACACTCGGGGGTAGTGCGCGCTGATGGCGA"
        snp1 = vcf_record.VcfRecord("ref\t3\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n")
        snp2 = vcf_record.VcfRecord("ref\t6\t.\tC\tG,T\t.\tPASS\tSVTYPE=SNP\n")
        snp2_1 = vcf_record.VcfRecord("ref\t6\t.\tC\tG\t.\tPASS\tSVTYPE=SNP\n")
        snp2_2 = vcf_record.VcfRecord("ref\t6\t.\tC\tT\t.\tPASS\tSVTYPE=SNP\n")
        indel1 = vcf_record.VcfRecord("ref\t5\t.\tTCG\tTC\t.\tPASS\tSVTYPE=INDEL\n")
        indel2 = vcf_record.VcfRecord("ref\t10\t.\tGG\tG\t.\tPASS\tSVTYPE=INDEL\n")
        snp3 = vcf_record.VcfRecord("ref\t20\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\n")
        vcf_records = [snp1, indel1, snp2, indel2, snp3]
        expected_gap4 = [snp1, indel1, snp2_1, snp2_2, snp3]
        expected_gap2 = [snp1, indel1, snp2_1, snp2_2, indel2, snp3]
        self.assertEqual(
            expected_gap2,
            vcf_clusterer.VcfClusterer._expand_alts_and_remove_duplicates_in_list(
                vcf_records, ref_seq, indel_gap=2
            ),
        )
        self.assertEqual(
            expected_gap4,
            vcf_clusterer.VcfClusterer._expand_alts_and_remove_duplicates_in_list(
                vcf_records, ref_seq, indel_gap=4
            ),
        )

    def test_cluster_vcf_record_list(self):
        """test _cluster_vcf_record_list"""
        record1 = vcf_record.VcfRecord(
            "ref_42\t11\tid_1\tA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        record2 = vcf_record.VcfRecord(
            "ref_42\t12\tid_2\tC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        record3 = vcf_record.VcfRecord(
            "ref_42\t15\tid_2\tC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        record4 = vcf_record.VcfRecord(
            "ref_42\t19\tid_2\tCCCCC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        record5 = vcf_record.VcfRecord(
            "ref_42\t23\tid_2\tC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        record_list = [record1, record2, record3, record4, record5]
        cluster1 = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record1, cluster_boundary_size=3
        )
        self.assertTrue(cluster1.add_vcf_record(record2))
        self.assertTrue(cluster1.add_vcf_record(record3))
        cluster2 = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record4, cluster_boundary_size=3
        )
        self.assertTrue(cluster2.add_vcf_record(record5))
        expected = [cluster1, cluster2]
        got = vcf_clusterer.VcfClusterer._cluster_vcf_record_list(
            record_list, cluster_boundary_size=3
        )
        self.assertEqual(expected, got)

        cluster1.cluster_boundary_size = 5
        self.assertTrue(cluster1.add_vcf_record(record4))
        self.assertTrue(cluster1.add_vcf_record(record5))
        got = vcf_clusterer.VcfClusterer._cluster_vcf_record_list(
            record_list, cluster_boundary_size=5
        )
        self.assertEqual([cluster1], got)

    def test_run_gramtools_merge(self):
        """test run using merge_method gramtools"""
        vcf_files = [
            os.path.join(data_dir, "run.gramtools_merge.1.vcf"),
            os.path.join(data_dir, "run.gramtools_merge.2.vcf"),
        ]
        ref_fasta = os.path.join(data_dir, "run.gramtools_merge.ref.fa")
        tmp_out = "tmp.vcf_clusterer.run.gramtools_merge.out.vcf"
        clusterer = vcf_clusterer.VcfClusterer(
            vcf_files,
            ref_fasta,
            tmp_out,
            source="source_name",
            max_alleles_per_cluster=8,
            cluster_boundary_size=1,
        )
        clusterer.run()
        expected_vcf = os.path.join(
            data_dir, "run.gramtools_merge.out.max_alleles_8.vcf"
        )
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)

        clusterer = vcf_clusterer.VcfClusterer(
            vcf_files,
            ref_fasta,
            tmp_out,
            source="source_name",
            max_alleles_per_cluster=100,
            cluster_boundary_size=1,
        )
        clusterer.run()
        expected_vcf = os.path.join(
            data_dir, "run.gramtools_merge.out.max_alleles_100.vcf"
        )
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_run_simple_merge(self):
        """test run using merge_method simple"""
        vcf_files = [
            os.path.join(data_dir, "run.simple_merge.1.vcf"),
            os.path.join(data_dir, "run.simple_merge.2.vcf"),
        ]
        ref_fasta = os.path.join(data_dir, "run.simple_merge.ref.fa")
        tmp_out = "tmp.vcf_clusterer.run.simple_merge.out.vcf"
        clusterer = vcf_clusterer.VcfClusterer(
            vcf_files, ref_fasta, tmp_out, source="source_name", merge_method="simple",
            cluster_boundary_size=1,
        )
        clusterer.run()
        expected_vcf = os.path.join(data_dir, "run.simple_merge.out.vcf")
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_run_gt_aware_merge(self):
        """test run using merge_method gt_aware"""
        vcf_files = [
            os.path.join(data_dir, "run.gt_aware_merge.1.vcf"),
            os.path.join(data_dir, "run.gt_aware_merge.2.vcf"),
        ]
        ref_fasta = os.path.join(data_dir, "run.gt_aware_merge.ref.fa")
        tmp_out = "tmp.vcf_clusterer.run.gt_aware_merge.out.vcf"
        clusterer = vcf_clusterer.VcfClusterer(
            vcf_files, ref_fasta, tmp_out, source="source_name", merge_method="gt_aware",
            cluster_boundary_size=1,
        )
        clusterer.run()
        expected_vcf = os.path.join(data_dir, "run.gt_aware_merge.out.vcf")
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)
