import copy
import os
import unittest

import pyfastaq

from cluster_vcf_records import vcf_record, vcf_record_cluster, vcf_clusterer

modules_dir = os.path.dirname(os.path.abspath(vcf_record_cluster.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "vcf_record_cluster")


class TestVcfRecordCluster(unittest.TestCase):
    def test_add_vcf_record_and_len(self):
        """test add_vcf_record and len"""
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

        cluster = vcf_record_cluster.VcfRecordCluster(cluster_boundary_size=3)
        self.assertEqual(0, len(cluster))
        self.assertTrue(cluster.add_vcf_record(record1))
        self.assertEqual(1, len(cluster))
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertEqual(2, len(cluster))
        self.assertTrue(cluster.add_vcf_record(record3))
        self.assertEqual(3, len(cluster))
        self.assertFalse(cluster.add_vcf_record(record4))
        self.assertEqual(3, len(cluster))
        cluster.cluster_boundary_size = 5
        self.assertTrue(cluster.add_vcf_record(record4))
        self.assertEqual(4, len(cluster))
        self.assertTrue(cluster.add_vcf_record(record5))
        self.assertEqual(5, len(cluster))

    def test_start_and_end(self):
        """test start_and_end"""
        cluster = vcf_record_cluster.VcfRecordCluster(cluster_boundary_size=3)
        self.assertEqual((None, None), cluster.start_and_end())
        record1 = vcf_record.VcfRecord(
            "ref_42\t11\tid_1\tA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        self.assertTrue(cluster.add_vcf_record(record1))
        self.assertEqual((10, 10), cluster.start_and_end())
        record2 = vcf_record.VcfRecord(
            "ref_42\t12\tid_2\tC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80"
        )
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertEqual((10, 11), cluster.start_and_end())

    def test_make_one_merged_vcf_record_for_gramtools(self):
        """test make_one_vcf_merged_record_for_gramtools"""
        #          12345678901234567890
        ref_seq = "AGCTATCTGCGTATTCGATC"
        record1 = vcf_record.VcfRecord(
            "ref\t10\t.\tC\tG\t42.42\t.\tSVTYPE=SNP\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record1, cluster_boundary_size=2
        )

        # If there's only one record, then we should just get that one back,
        # but should have FILTER = PASS
        expected = copy.copy(record1)
        expected.FILTER = {"PASS"}
        self.assertEqual(
            expected, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

        # Test one record, but two alleles
        record2 = vcf_record.VcfRecord(
            "ref\t12\t.\tT\tC,A\t42.42\tPASS\tSVTYPE=SNP\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record2, cluster_boundary_size=2
        )
        self.assertEqual(
            record2, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

        # Two records, SNPs at different positions
        self.assertTrue(cluster.add_vcf_record(record1))
        expected = vcf_record.VcfRecord(
            "ref\t10\t.\tCGT\tCGA,CGC,GGA,GGC,GGT\t.\tPASS\tSVTYPE=COMPLEX"
        )
        self.assertEqual(
            expected, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

        # Add in a deletion to the previous two SNPs, starting a earlier and ending
        # later, to also check flanking regions get added.
        record3 = vcf_record.VcfRecord(
            "ref\t8\t.\tTGCGTAT\tT\t42.42\tPASS\tSVTYPE=COMPLEX\tGT\t1/1"
        )
        self.assertTrue(cluster.add_vcf_record(record3))
        expected = vcf_record.VcfRecord(
            "ref\t8\t.\tTGCGTAT\tT,TGCGAAT,TGCGCAT,TGGGAAT,TGGGCAT,TGGGTAT\t.\tPASS\tSVTYPE=COMPLEX"
        )
        self.assertEqual(
            expected, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

        # Add an indel to all the previous variants.
        record4 = vcf_record.VcfRecord(
            "ref\t9\t.\tGCGT\tGAAAAC\t42.42\tPASS\tSVTYPE=COMPLEX\tGT\t1/1"
        )
        self.assertTrue(cluster.add_vcf_record(record4))
        expected = vcf_record.VcfRecord(
            "ref\t8\t.\tTGCGTAT\tT,TGAAAACAT,TGCGAAT,TGCGCAT,TGGGAAT,TGGGCAT,TGGGTAT\t.\tPASS\tSVTYPE=COMPLEX"
        )
        self.assertEqual(
            expected, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

        # Test insertion next to SNP
        record5 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCG\t42.42\tPASS\tSVTPYPE=INDEL\tGT\t1/1"
        )
        record6 = vcf_record.VcfRecord(
            "ref\t4\t.\tT\tA\t42.42\tPASS\tSVTPYPE=SNP\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record5, cluster_boundary_size=1
        )
        self.assertTrue(cluster.add_vcf_record(record6))
        expected = vcf_record.VcfRecord(
            "ref\t3\t.\tCT\tCA,CGA,CGT\t.\tPASS\tSVTYPE=COMPLEX"
        )
        self.assertEqual(
            expected, cluster.make_one_merged_vcf_record_for_gramtools(ref_seq)
        )

    def test_SeveralComplexRecords_CorrectClustering(self):
        """
        Here's an example where we have two complex records followed by a SNP.
        This unit test i)tests cluster chaining; ii)illustrates importance of record POS ordering.

        Assertions 1 and 2:
            Combining the first and third records yields two distinct clusters;
            combining them all, in order, together gives a single cluster.

            This is because extending the first with the second (NB: in that order!!) extends the cluster boundaries
            to make the third record overlap.

        Assertion 3:
            Here I show the importance of order. First note that we only look at the last cluster in the list of clusters
            to determine if a record gets added to that, or makes a new cluster.

            Adding record 2 AFTER records 1 and 3, the fact that record 2 overlaps with record 1 is lost;
            and so is the bridging of all three.
        """
        # ref_seq = 'TTCGGCAGGCATT' # For the record.
        record_1 = vcf_record.VcfRecord(
            "ref\t3\t.\tCGGCA\tC\t.\tPASS\tSVTPYPE=COMPLEX\tGT\t1/1"
        )
        record_2 = vcf_record.VcfRecord(
            "ref\t5\t.\tGCAGGCA\tGCTGGCT\t.\tPASS\tSVTPYPE=COMPLEX\tGT\t1/1"
        )
        record_3 = vcf_record.VcfRecord(
            "ref\t10\t.\tC\tA\t42.42\tPASS\tSVTPYPE=SNP\tGT\t1/1"
        )

        Two_clusters = vcf_clusterer.VcfClusterer._cluster_vcf_record_list(
            [record_1, record_3], cluster_boundary_size=1
        )
        self.assertEqual(len(Two_clusters), 2)

        One_cluster = vcf_clusterer.VcfClusterer._cluster_vcf_record_list(
            [record_1, record_2, record_3], cluster_boundary_size=1
        )

        self.assertEqual(len(One_cluster), 1)

        Three_clusters = vcf_clusterer.VcfClusterer._cluster_vcf_record_list(
            [record_1, record_3, record_2], cluster_boundary_size=1
        )
        self.assertEqual(len(Three_clusters), 3)

    def test_make_simple_merged_vcf_with_no_combinations(self):
        """test make_simple_merged_vcf_with_no_combinations"""
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTAGGTCAG")
        record1 = vcf_record.VcfRecord(
            "ref\t2\t.\tG\tAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCAT\t21.4018\t.\tINDEL;IDV=2;IMF=0.0338983;DP=59;VDB=0.18;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60\tGT:PL\t1/1:48,6,0"
        )
        record3 = vcf_record.VcfRecord(
            "ref\t8\t.\tT\tA\t21.4018\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(vcf_record=record1, cluster_boundary_size=8)
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertTrue(cluster.add_vcf_record(record3))
        cluster.make_simple_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t2\t.\tGCTAGGT\tAACATTAGGA\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        self.assertEqual(expected, cluster[0])

        # add in record 1 again, then try merging. The merging should not work
        self.assertTrue(cluster.add_vcf_record(record1))
        cluster.make_simple_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(2, len(cluster))
        self.assertEqual(expected, cluster[0])
        self.assertEqual(record1, cluster[1])

        # Test insertion next to SNP
        # (same example as in test_make_one_merged_vcf_record_for_gramtools,
        # byt now we don't get combinations, just the two ALT changes from
        # the VCF records)
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTATCTGCGTATTCGATC")
        record1 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCG\t42.42\tPASS\tSVTPYPE=INDEL\tGT\t1/1"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t4\t.\tT\tA\t42.42\tPASS\tSVTPYPE=SNP\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record1, cluster_boundary_size=1
        )
        self.assertTrue(cluster.add_vcf_record(record2))
        cluster.make_simple_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t3\t.\tCT\tCGA\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        self.assertEqual(expected, cluster[0])

    def test_make_simple_gt_aware_merged_vcf_with_no_combinations(self):
        """test make_simple_merged_vcf_with_no_combinations"""
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTAGGTCAG")
        record1 = vcf_record.VcfRecord(
            "ref\t2\t.\tG\tAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCAT\t21.4018\t.\tINDEL;IDV=2;IMF=0.0338983;DP=59;VDB=0.18;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60\tGT:PL\t0/0:48,6,0"
        )
        record3 = vcf_record.VcfRecord(
            "ref\t8\t.\tT\tA\t21.4018\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(vcf_record=record1, cluster_boundary_size=8)
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertTrue(cluster.add_vcf_record(record3))
        cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t2\t.\tGCTAGGT\tAACTAGGA\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        self.assertEqual(expected, cluster[0])

        # add in record 1 again, then try merging. The merging should not work
        self.assertTrue(cluster.add_vcf_record(record1))
        cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(2, len(cluster))
        self.assertEqual(expected, cluster[0])
        self.assertEqual(record1, cluster[1])

        # repeat with different combinations of ref/alt calls
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTAGGTCAG")
        record1 = vcf_record.VcfRecord(
            "ref\t2\t.\tG\tAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t0/0:255,163,0"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCAT\t21.4018\t.\tINDEL;IDV=2;IMF=0.0338983;DP=59;VDB=0.18;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60\tGT:PL\t1/1:48,6,0"
        )
        record3 = vcf_record.VcfRecord(
            "ref\t8\t.\tT\tA\t21.4018\t.\tSVTYPE=MERGED\tGT\t0/0"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(vcf_record=record1, cluster_boundary_size=8)
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertTrue(cluster.add_vcf_record(record3))
        cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t2\t.\tGCTAGGT\tGCATTAGGT\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        self.assertEqual(expected, cluster[0])

        # test where all calls are ref
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTAGGTCAG")
        record1 = vcf_record.VcfRecord(
            "ref\t2\t.\tG\tAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t0/0:255,163,0"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCAT\t21.4018\t.\tINDEL;IDV=2;IMF=0.0338983;DP=59;VDB=0.18;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60\tGT:PL\t0/0:48,6,0"
        )
        record3 = vcf_record.VcfRecord(
            "ref\t8\t.\tT\tA\t21.4018\t.\tSVTYPE=MERGED\tGT\t0/0"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(vcf_record=record1, cluster_boundary_size=8)
        self.assertTrue(cluster.add_vcf_record(record2))
        self.assertTrue(cluster.add_vcf_record(record3))
        cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t2\t.\tGCTAGGT\tAACATTAGGA\t.\t.\tSVTYPE=MERGED\tGT\t0/0"
        )
        self.assertEqual(expected, cluster[0])

        # Test insertion next to SNP
        # (same example as in test_make_one_merged_vcf_record_for_gramtools,
        # byt now we don't get combinations, just the two ALT changes from
        # the VCF records)
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTATCTGCGTATTCGATC")
        record1 = vcf_record.VcfRecord(
            "ref\t3\t.\tC\tCG\t42.42\tPASS\tSVTPYPE=INDEL\tGT\t1/1"
        )
        record2 = vcf_record.VcfRecord(
            "ref\t4\t.\tT\tA\t42.42\tPASS\tSVTPYPE=SNP\tGT\t1/1"
        )
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=record1, cluster_boundary_size=1
        )
        self.assertTrue(cluster.add_vcf_record(record2))
        cluster.make_simple_gt_aware_merged_vcf_with_no_combinations(ref_seq)
        self.assertEqual(1, len(cluster))
        expected = vcf_record.VcfRecord(
            "ref\t3\t.\tCT\tCGA\t.\t.\tSVTYPE=MERGED\tGT\t1/1"
        )
        self.assertEqual(expected, cluster[0])

    def test_make_separate_indels_and_one_alt_with_all_snps_no_combinations(self):
        """test make_separate_indels_and_one_alt_with_all_snps_no_combinations"""
        ref_seq = pyfastaq.sequences.Fasta("ref", "AGCTAGGTCAG")
        snp1 = vcf_record.VcfRecord("ref\t4\t.\tT\tA\t.\t.\t.\t.")
        snp2 = vcf_record.VcfRecord("ref\t9\t.\tC\tG\t.\t.\t.\t.")
        complex_var = vcf_record.VcfRecord(
            "ref\t3\t.\tCTAGGTCA\tCTATTGGTCA\t.\t.\t.\t."
        )
        deletion = vcf_record.VcfRecord("ref\t3\t.\tCTAGGTCA\tG\t.\t.\t.\t.")
        insertion = vcf_record.VcfRecord("ref\t5\t.\tA\tATT\t.\t.\t.\t.")
        cluster = vcf_record_cluster.VcfRecordCluster(
            vcf_record=deletion, cluster_boundary_size=1
        )
        self.assertTrue(cluster.add_vcf_record(snp1))
        self.assertTrue(cluster.add_vcf_record(snp2))
        self.assertTrue(cluster.add_vcf_record(insertion))
        self.assertTrue(cluster.add_vcf_record(complex_var))
        got = cluster.make_separate_indels_and_one_alt_with_all_snps_no_combinations(
            ref_seq
        )
        expected = vcf_record.VcfRecord(
            "ref\t3\t.\tCTAGGTCA\tCAAGGTGA,CTATTGGTCA,G\t.\tPASS\t."
        )
        self.assertEqual(expected, got)
