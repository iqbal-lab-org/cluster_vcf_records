import itertools
import unittest

import pyfastaq

from cluster_vcf_records import vcf_record

class TestVcfRecord(unittest.TestCase):
    def test_VcfRecord_constructor(self):
        '''test VcfRecord constructor'''
        line = 'ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(record.CHROM, 'ref_42')
        self.assertEqual(record.POS, 10)
        self.assertEqual(record.ID, 'id_foo')
        self.assertEqual(record.REF, 'A')
        self.assertEqual(record.ALT, ['G'])
        self.assertEqual(record.QUAL, 42.42)
        self.assertEqual(record.FILTER, {'PASS'})
        self.assertEqual(record.INFO, {'KMER': '31', 'SVLEN': '0', 'SVTYPE': 'SNP'})
        self.assertEqual(record.FORMAT, {'GT': '1/1', 'COV': '0,52', 'GT_CONF': '39.80'})

        line = 'ref_42\t11\tid_foo\tA\tG,TC\t.\tFilter1;Filter2\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(record.QUAL, None)
        self.assertEqual(record.ALT, ['G', 'TC'])
        self.assertEqual(record.FILTER, {'Filter1', 'Filter2'})

        line = 'ref_42\t11\tid_foo\tA\tG,TC\t.\t.\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(record.FILTER, set())

    def test_str(self):
        '''test __str__'''
        line = 'ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(line, str(record))

        line = 'ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\t.'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(line, str(record))

        line = 'ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\tFOO;KMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(line, str(record))


    def test_ref_string_matches_ref_sequence(self):
        '''test ref_string_matches_ref_sequence'''
        record = vcf_record.VcfRecord('ref_name\t-1\t.\tAGT\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.ref_string_matches_ref_sequence('AG'))
        record = vcf_record.VcfRecord('ref_name\t3\t.\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertTrue(record.ref_string_matches_ref_sequence('GCATG'))
        self.assertFalse(record.ref_string_matches_ref_sequence('GCxTG'))
        record = vcf_record.VcfRecord('ref_name\t3\t.\tAGT\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertTrue(record.ref_string_matches_ref_sequence('GCAGT'))
        self.assertFalse(record.ref_string_matches_ref_sequence('GCAGC'))
        self.assertFalse(record.ref_string_matches_ref_sequence('GCAG'))
        self.assertFalse(record.ref_string_matches_ref_sequence('GCA'))
        self.assertFalse(record.ref_string_matches_ref_sequence('GA'))


    def test_ref_string_matches_dict_of_ref_sequences(self):
        '''test ref_string_matches_dict_of_ref_sequences'''
        record = vcf_record.VcfRecord('ref1\t3\t.\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        ref_seqs = {'ref1': 'GTACG', 'ref2': 'TTTTT'}
        self.assertTrue(record.ref_string_matches_dict_of_ref_sequences(ref_seqs))
        record = vcf_record.VcfRecord('ref2\t3\t.\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.ref_string_matches_dict_of_ref_sequences(ref_seqs))
        record = vcf_record.VcfRecord('ref3\t3\t.\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.ref_string_matches_dict_of_ref_sequences(ref_seqs))


    def test_remove_asterisk_alts(self):
        '''test remove_asterisk_alts'''
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        record.remove_asterisk_alts()
        self.assertEqual(['G'], record.ALT)
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tG,*\tPASS\tSVTYPE=SNP\tGT\t1/1')
        record.remove_asterisk_alts()
        self.assertEqual(['G'], record.ALT)
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\t*\tPASS\tSVTYPE=SNP\tGT\t1/1')
        record.remove_asterisk_alts()
        self.assertEqual([], record.ALT)


    def test_is_snp(self):
        '''test is_snp'''
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertTrue(record.is_snp())
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tG,T\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertTrue(record.is_snp())
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tGT\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.is_snp())
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tA\tG,GT\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.is_snp())
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tTG\tG\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.is_snp())
        record = vcf_record.VcfRecord('ref.3\t8\tid5\tTG\t.\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertFalse(record.is_snp())


    def test_is_homozygous(self):
        '''test is_homozygous'''
        record = vcf_record.VcfRecord('ref\t1\tid1\tA\tG\t42\tPASS\tSVTYPE=SNP\tGT\t0/0')
        self.assertTrue(record.is_homozygous())
        record = vcf_record.VcfRecord('ref\t1\tid1\tA\tG\t42\tPASS\tSVTYPE=SNP\tGT\t1/1')
        self.assertTrue(record.is_homozygous())
        record = vcf_record.VcfRecord('ref\t1\tid1\tA\tG,T\t42\tPASS\tSVTYPE=SNP\tGT\t2/2')
        self.assertTrue(record.is_homozygous())
        record = vcf_record.VcfRecord('ref\t1\tid1\tA\tG,T\t42\tPASS\tSVTYPE=SNP\tGT\t1/2')
        self.assertFalse(record.is_homozygous())
        record = vcf_record.VcfRecord('ref\t1\tid1\tA\tG,T\t42\tPASS\tSVTYPE=SNP\tGT\t./.')
        self.assertFalse(record.is_homozygous())


    def test_set_format_key_value(self):
        '''test set_format_key_value'''
        record = vcf_record.VcfRecord('ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\tFOO;KMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80')
        self.assertFalse('FOO' in record.FORMAT)
        record.set_format_key_value('FOO', 'bar')
        self.assertTrue('FOO' in record.FORMAT)
        self.assertEqual('bar', record.FORMAT['FOO'])
        record.set_format_key_value('FOO', 'baz')
        self.assertEqual('baz', record.FORMAT['FOO'])


    def test_intersects(self):
        '''test intersects'''
        record1 = vcf_record.VcfRecord('ref_42\t11\t.\tA\tG\t42.42\tPASS\t.')
        record2 = vcf_record.VcfRecord('ref_42\t12\t.\tC\tT\t42.42\tPASS\t.')
        record3 = vcf_record.VcfRecord('ref_43\t12\t.\tC\tT\t42.42\tPASS\t.')
        record4 = vcf_record.VcfRecord('ref_42\t11\t.\tCT\tT\t42.42\tPASS\t.')
        self.assertTrue(record1.intersects(record1))
        self.assertTrue(record2.intersects(record2))
        self.assertFalse(record1.intersects(record2))
        self.assertFalse(record2.intersects(record1))
        self.assertFalse(record3.intersects(record2))
        self.assertFalse(record2.intersects(record3))
        self.assertTrue(record2.intersects(record4))
        self.assertTrue(record4.intersects(record2))


    def test_merge(self):
        '''test merge'''
        ref_seq = pyfastaq.sequences.Fasta('ref', 'AGCTAGGTCAG')
        record1 = vcf_record.VcfRecord('wrong_ref\t3\t.\tC\tA\t228\t.\t.\t.')
        record2 = vcf_record.VcfRecord('ref\t1\t.\tAG\tGAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record3 = vcf_record.VcfRecord('ref\t2\t.\tG\tAA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record4 = vcf_record.VcfRecord('ref\t3\t.\tC\tCAT\t21.4018\t.\tINDEL;IDV=2;IMF=0.0338983;DP=59;VDB=0.18;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60\tGT:PL\t1/1:48,6,0')
        record5 = vcf_record.VcfRecord('ref\t7\t.\tG\tC\t21.4018\t.\t.\t.\t.')

        self.assertIsNone(record1.merge(record2, ref_seq))
        self.assertIsNone(record2.merge(record1, ref_seq))
        self.assertIsNone(record1.merge(record3, ref_seq))
        self.assertIsNone(record3.merge(record1, ref_seq))

        got = record3.merge(record4, ref_seq)
        expected = vcf_record.VcfRecord('ref\t2\t.\tGC\tAACAT\t.\t.\tSVTYPE=MERGED\tGT\t1/1')
        self.assertEqual(expected, got)
        got = record4.merge(record3, ref_seq)
        expected = vcf_record.VcfRecord('ref\t2\t.\tGC\tAACAT\t.\t.\tSVTYPE=MERGED\tGT\t1/1')
        self.assertEqual(expected, got)

        got = record4.merge(record5, ref_seq)
        expected = vcf_record.VcfRecord('ref\t3\t.\tCTAGG\tCATTAGC\t.\t.\tSVTYPE=MERGED\tGT\t1/1')
        self.assertEqual(expected, got)


    def test_add_flanking_seqs(self):
        '''test add_flanking_seqs'''
        #                                          01234567890
        ref_seq = pyfastaq.sequences.Fasta('ref', 'AGCTAGGTCAG')
        record = vcf_record.VcfRecord('ref\t4\t.\tT\tG,TC\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record.add_flanking_seqs(ref_seq, 1, 6)
        self.assertEqual(1, record.POS)
        self.assertEqual(6, record.ref_end_pos())
        self.assertEqual(['GCGAGG', 'GCTCAGG'], record.ALT)
        with self.assertRaises(vcf_record.Error):
            record.add_flanking_seqs(ref_seq, 2, 10)
        with self.assertRaises(vcf_record.Error):
            record.add_flanking_seqs(ref_seq, 0, 5)


    def test_merge_by_adding_new_alts(self):
        '''test merge_by_adding_new_alts'''
        #                          ref coords:     01234567890
        ref_seq = pyfastaq.sequences.Fasta('ref', 'AGTGAGGTCAG')
        record1 = vcf_record.VcfRecord('ref\t3\t.\tT\tG\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record2 = vcf_record.VcfRecord('ref\t3\t.\tT\tG\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record1.merge_by_adding_new_alts(record2, ref_seq)
        self.assertEqual((2, 'T', ['G']), (record1.POS, record1.REF, record1.ALT))

        record2 = vcf_record.VcfRecord('ref\t3\t.\tT\tA\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record1.merge_by_adding_new_alts(record2, ref_seq)
        self.assertEqual((2, 'T', ['G', 'A']), (record1.POS, record1.REF, record1.ALT))
        record1.merge_by_adding_new_alts(record2, ref_seq)
        self.assertEqual((2, 'T', ['G', 'A']), (record1.POS, record1.REF, record1.ALT))

        record2 = vcf_record.VcfRecord('ref\t4\t.\tG\tT,C\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record1.merge_by_adding_new_alts(record2, ref_seq)
        self.assertEqual((2, 'TG', ['GG', 'AG', 'TT', 'TC']), (record1.POS, record1.REF, record1.ALT))

        record1 = vcf_record.VcfRecord('ref\t7\t.\tGT\tG,GAC\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record2 = vcf_record.VcfRecord('ref\t3\t.\tT\tG\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record1.merge_by_adding_new_alts(record2, ref_seq)
        self.assertEqual((2, 'TGAGGT', ['TGAGG', 'TGAGGAC', 'GGAGGT']), (record1.POS, record1.REF, record1.ALT))


    def test_ref_end_pos(self):
        '''test ref_end_pos'''
        line = 'ref_42\t11\tid_foo\tA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(10, record.ref_end_pos())

        line = 'ref_42\t11\tid_foo\tAA\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(11, record.ref_end_pos())

        line = 'ref_42\t11\tid_foo\tAAG\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(12, record.ref_end_pos())

        line = 'ref_42\t11\tid_foo\t.\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n'
        record = vcf_record.VcfRecord(line)
        self.assertEqual(10, record.ref_end_pos())


    def test_remove_useless_start_nucleotides(self):
        '''test remove_useless_start_nucleotides'''
        record = vcf_record.VcfRecord('ref\t9002\t.\tGT\tGTT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        expected = vcf_record.VcfRecord('ref\t9003\t.\tT\tTT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record.remove_useless_start_nucleotides()
        self.assertEqual(expected, record)

        record = vcf_record.VcfRecord('ref\t9004\t.\tGTTT\tGT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        expected = vcf_record.VcfRecord('ref\t9005\t.\tTTT\tT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record.remove_useless_start_nucleotides()
        self.assertEqual(expected, record)

        record = vcf_record.VcfRecord('ref\t9006\t.\tG\tT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        expected = vcf_record.VcfRecord('ref\t9006\t.\tG\tT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record.remove_useless_start_nucleotides()
        self.assertEqual(expected, record)

        record = vcf_record.VcfRecord('ref\t9007\t.\t.\tT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        expected = vcf_record.VcfRecord('ref\t9007\t.\t.\tT\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0')
        record.remove_useless_start_nucleotides()
        self.assertEqual(expected, record)


    def test_near_to_position(self):
        '''test near_to_position'''
        record = vcf_record.VcfRecord('ref_42\t11\tid_foo\tAAC\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertFalse(record.near_to_position(8, 1))
        self.assertTrue(record.near_to_position(8, 2))
        self.assertTrue(record.near_to_position(9, 1))
        self.assertTrue(record.near_to_position(10, 1))
        self.assertTrue(record.near_to_position(11, 1))
        self.assertTrue(record.near_to_position(12, 1))
        self.assertTrue(record.near_to_position(13, 1))
        self.assertFalse(record.near_to_position(14, 1))
        self.assertTrue(record.near_to_position(14, 2))


    def test_inferred_var_seqs_plus_flanks(self):
        '''test inferred_var_seqs_plus_flanks'''
        #          12345678901234567890123456
        ref_seq = 'CGCGCGAGCGTGTAGAGTGCCAGACT'
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((2, ['C', 'A']), record.inferred_var_seqs_plus_flanks(ref_seq, 0))
        self.assertEqual((1, ['GCG', 'GAG']), record.inferred_var_seqs_plus_flanks(ref_seq, 1))
        self.assertEqual((0, ['CGCGC', 'CGAGC']), record.inferred_var_seqs_plus_flanks(ref_seq, 2))
        self.assertEqual((0, ['CGCGCG', 'CGAGCG']), record.inferred_var_seqs_plus_flanks(ref_seq, 3))

        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA,GT\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((0, ['CGCGCG', 'CGAGCG', 'CGGTGCG']), record.inferred_var_seqs_plus_flanks(ref_seq, 3))

        record = vcf_record.VcfRecord('ref\t24\tid_foo\tA\tT\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((23, ['A', 'T']), record.inferred_var_seqs_plus_flanks(ref_seq, 0))
        self.assertEqual((22, ['GAC', 'GTC']), record.inferred_var_seqs_plus_flanks(ref_seq, 1))
        self.assertEqual((21, ['AGACT', 'AGTCT']), record.inferred_var_seqs_plus_flanks(ref_seq, 2))
        self.assertEqual((20, ['CAGACT', 'CAGTCT']), record.inferred_var_seqs_plus_flanks(ref_seq, 3))
        self.assertEqual((19, ['CCAGACT', 'CCAGTCT']), record.inferred_var_seqs_plus_flanks(ref_seq, 4))

        record = vcf_record.VcfRecord('ref\t10\tid_foo\tGT\tG\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((9, ['GT', 'G']), record.inferred_var_seqs_plus_flanks(ref_seq, 0))
        self.assertEqual((8, ['CGTG', 'CGG']), record.inferred_var_seqs_plus_flanks(ref_seq, 1))
        self.assertEqual((7, ['GCGTGT', 'GCGGT']), record.inferred_var_seqs_plus_flanks(ref_seq, 2))
        self.assertEqual((6, ['AGCGTGTA', 'AGCGGTA']), record.inferred_var_seqs_plus_flanks(ref_seq, 3))

        record = vcf_record.VcfRecord('ref\t10\tid_foo\tG\tGAA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((9, ['G', 'GAA']), record.inferred_var_seqs_plus_flanks(ref_seq, 0))
        self.assertEqual((8, ['CGT', 'CGAAT']), record.inferred_var_seqs_plus_flanks(ref_seq, 1))
        self.assertEqual((7, ['GCGTG', 'GCGAATG']), record.inferred_var_seqs_plus_flanks(ref_seq, 2))

        record = vcf_record.VcfRecord('ref\t10\tid_foo\tGTGT\tATA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:0,52:39.80\n')
        self.assertEqual((9, ['GTGT', 'ATA']), record.inferred_var_seqs_plus_flanks(ref_seq, 0))
        self.assertEqual((8, ['CGTGTA', 'CATAA']), record.inferred_var_seqs_plus_flanks(ref_seq, 1))
        self.assertEqual((7, ['GCGTGTAG', 'GCATAAG']), record.inferred_var_seqs_plus_flanks(ref_seq, 2))
        self.assertEqual((6, ['AGCGTGTAGA', 'AGCATAAGA']), record.inferred_var_seqs_plus_flanks(ref_seq, 3))


    def test_total_coverage(self):
        '''test total_coverage'''
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:GT_CONF\t1/1:39.80\n')
        self.assertEqual(None, record.total_coverage())
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:COV:GT_CONF\t1/1:1,2,39:39.80\n')
        self.assertEqual(42, record.total_coverage())


    def test_called_alts_from_genotype(self):
        '''test called_alts_from_genotype'''
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT_CONF\t39.80\n')
        self.assertEqual(None, record.called_alts_from_genotype())
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA,G,T\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:GT_CONF\t0/0:39.80\n')
        self.assertEqual({'C'}, record.called_alts_from_genotype())
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA,G,T\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:GT_CONF\t1/2:39.80\n')
        self.assertEqual({'A', 'G'}, record.called_alts_from_genotype())
        record = vcf_record.VcfRecord('ref\t3\tid_foo\tC\tA,G,T\t42.42\tPASS\tKMER=31;SVLEN=0;SVTYPE=SNP\tGT:GT_CONF\t1/1:39.80\n')
        self.assertEqual({'A'}, record.called_alts_from_genotype())


    def test_is_the_same_indel(self):
        '''test is_the_same_indel'''
        #          123456789012345678901234567890
        ref_seq = 'ACACTCGGGTAGTCAGCGCGCTGATGGCGA'
        snp = vcf_record.VcfRecord('ref\t3\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n')

        # test homopolymer runs
        deletions = [
            vcf_record.VcfRecord('ref\t5\t.\tTCGGG\tTCGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCGGG\tCGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCGG\tCG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCG\tC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t7\t.\tGGG\tGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t7\t.\tGG\tG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t8\t.\tGG\tG\t.\tPASS\tSVTYPE=SNP\n'),
        ]
        insertions = [
            vcf_record.VcfRecord('ref\t5\t.\tTC\tTCG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tC\tCG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCG\tCGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCGG\tCGGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t6\t.\tCGGG\tCGGGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t7\t.\tG\tGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t8\t.\tG\tGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t8\t.\tGG\tGGG\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t9\t.\tG\tGG\t.\tPASS\tSVTYPE=SNP\n'),
        ]

        del1 = vcf_record.VcfRecord('ref\t6\t.\tCGGG\tCG\t.\tPASS\tSVTYPE=SNP\n')
        ins1 = vcf_record.VcfRecord('ref\t5\t.\tTC\tTCA\t.\tPASS\tSVTYPE=SNP\n')


        for vcf1, vcf2 in itertools.permutations(deletions, 2):
            self.assertTrue(vcf1.is_the_same_indel(vcf2, ref_seq))

        for vcf1, vcf2 in itertools.permutations(insertions, 2):
            self.assertTrue(vcf1.is_the_same_indel(vcf2, ref_seq))

        for vcf in deletions:
            self.assertFalse(vcf.is_the_same_indel(del1, ref_seq))
            self.assertFalse(vcf.is_the_same_indel(snp, ref_seq))

        for vcf in insertions:
            self.assertFalse(vcf.is_the_same_indel(ins1, ref_seq))
            self.assertFalse(vcf.is_the_same_indel(snp, ref_seq))

        # test dinucleotide repeats
        deletions = [
            vcf_record.VcfRecord('ref\t15\t.\tAGCGCGC\tAGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t15\t.\tAGCGC\tAGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t15\t.\tAGC\tA\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t17\t.\tCGCGC\tCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t17\t.\tCGC\tC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t19\t.\tCGC\tC\t.\tPASS\tSVTYPE=SNP\n'),
        ]

        insertions = [
            vcf_record.VcfRecord('ref\t15\t.\tA\tAGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t15\t.\tAGC\tAGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t15\t.\tAGCGC\tAGCGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t15\t.\tAGCGCGC\tAGCGCGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t17\t.\tC\tCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t17\t.\tCGC\tCGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t17\t.\tCGCGC\tCGCGCGC\t.\tPASS\tSVTYPE=SNP\n'),
            vcf_record.VcfRecord('ref\t19\t.\tCGC\tCGCGC\t.\tPASS\tSVTYPE=SNP\n'),
        ]

        for vcf1, vcf2 in itertools.permutations(deletions, 2):
            self.assertTrue(vcf1.is_the_same_indel(vcf2, ref_seq))

        for vcf1, vcf2 in itertools.permutations(insertions, 2):
            self.assertTrue(vcf1.is_the_same_indel(vcf2, ref_seq))


    def test_to_record_per_alt(self):
        '''test to_record_per_alt'''
        vcf = vcf_record.VcfRecord('ref\t42\t.\tA\tC\t.\tPASS\tSVTYPE=SNP\n')
        self.assertEqual([vcf], vcf.to_record_per_alt())
        vcf = vcf_record.VcfRecord('ref\t42\t.\tA\tC,TC\t.\tPASS\tSVTYPE=SNP\n')
        vcf_c = vcf_record.VcfRecord('ref\t42\t.\tA\tC\t.\tPASS\tSVTYPE=SNP\n')
        vcf_tc = vcf_record.VcfRecord('ref\t42\t.\tA\tTC\t.\tPASS\tSVTYPE=SNP\n')
        self.assertEqual([vcf_c, vcf_tc], vcf.to_record_per_alt())


    def test_split_into_snps(self):
        '''test split_into_snps'''
        vcf = vcf_record.VcfRecord('ref\t15\t.\tA\tAG\t.\tPASS\tSVTYPE=SNP\n')
        expect = [vcf]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tA\tAG\t.\tPASS\tSVTYPE=SNP\n')
        expect = [vcf]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n')
        expect = [vcf]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tAT\tGC\t.\tPASS\tSVTYPE=SNP\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n')
        snp2 = vcf_record.VcfRecord('ref\t16\t.\tT\tC\t.\tPASS\tSVTYPE=SNP\n')
        expect = [snp1, snp2]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tAT\tGC\t.\tPASS\tSVTYPE=SNP\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n')
        snp2 = vcf_record.VcfRecord('ref\t16\t.\tT\tC\t.\tPASS\tSVTYPE=SNP\n')
        expect = [snp1, snp2]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGC,ATAA,GTGA,ATGG\t.\tPASS\tSVTYPE=SNP\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGG,ATAA,GTGG,ATGG\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGG,ATAA,GTGG,ATGG\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\tGT\t0/2\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGG,ATAA,GTGG,ATGG\t.\tPASS\tSVTYPE=SNP\tGT\t0/2\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\tGT\t0/0\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGG,ATAA,GTGG,ATGG\t.\tPASS\tSVTYPE=SNP\tGT\t1/2\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\tGT\t1/2\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)

        vcf = vcf_record.VcfRecord('ref\t15\t.\tATGC\tGTGA,ATAA,GTGG,ATGG\t.\tPASS\tSVTYPE=SNP\tGT\t1/2\n')
        snp1 = vcf_record.VcfRecord('ref\t15\t.\tA\tG\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp2 = vcf_record.VcfRecord('ref\t17\t.\tG\tA\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        snp3 = vcf_record.VcfRecord('ref\t18\t.\tC\tA,G\t.\tPASS\tSVTYPE=SNP\tGT\t0/1\n')
        expect = [snp1, snp2, snp3]
        got = vcf.split_into_snps()
        self.assertEqual(expect, got)
