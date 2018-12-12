import copy
import operator

class Error (Exception): pass


class VcfRecord:
    def __init__(self, line):
        '''Constructs VcfRecord from a line of a VCF file.
        Assumes only one sample in the file'''
        assert not line.startswith('#')
        fields = line.rstrip().split('\t')
        # #CHROM          POS     ID                      REF     ALT     QUAL    FILTER  INFO                            FORMAT          2.2.2.1
        # NC_000962.3     1977    UNION_BC_k31_var_120    A       G       .       PASS    KMER=31;SVLEN=0;SVTYPE=SNP      GT:COV:GT_CONF  1/1:0,52:39.80
        try:
            self.CHROM = fields[0]
            self.POS = int(fields[1]) - 1
            self.ID = fields[2]
            self.REF = fields[3]
            self.ALT = fields[4].split(',')
            self.QUAL = fields[5]
            self.FILTER = fields[6]
            INFO = fields[7]
        except:
            raise Error('Error reading line of vcf file:' + line)

        try:
            self.QUAL = float(self.QUAL)
        except:
            self.QUAL = None

        self.INFO = {}
        if INFO != '.':
            info_fields = INFO.split(';')
            for field in info_fields:
                if '=' in field:
                    key, value = field.split('=')
                    self.INFO[key] = value
                else:
                    self.INFO[field] = None

        if len(fields) == 10:
            FORMAT = fields[8]
            FORMAT_VALS = fields[9]
            self.format_keys = FORMAT.split(':')
            format_vals = FORMAT_VALS.split(':')
            self.FORMAT = dict(zip(self.format_keys, format_vals))
        else:
            self.format_keys = self.FORMAT = None


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __repr__(self):
        if len(self.INFO) == 0:
            info_string = '.'
        else:
            info_fields = []
            for x in sorted(self.INFO):
                if self.INFO[x] is None:
                    info_fields.append(x)
                else:
                    info_fields.append(x + '=' + self.INFO[x])
            info_string = ';'.join(info_fields)

        fields = [
            self.CHROM,
            str(self.POS + 1),
            self.ID,
            self.REF,
            ','.join(self.ALT),
            '.' if self.QUAL is None else str(self.QUAL),
            self.FILTER,
            info_string,
        ]

        if self.format_keys is not None:
            format_string = ':'.join(self.format_keys)
            format_values = ':'.join([self.FORMAT[x] for x in self.format_keys])
            fields.extend([format_string, format_values])

        return '\t'.join(fields)


    def remove_asterisk_alts(self):
        '''Removes "*" from alts, if there. Warning:
        this could mean that the records ends up with NO alts!'''
        self.ALT = [x for x in self.ALT if x != '*']


    def is_snp(self):
        '''Returns true iff this variant is a SNP'''
        nucleotides = {'A', 'C', 'G', 'T'}
        return len(self.REF) == 1 and self.REF in nucleotides and set(self.ALT).issubset(nucleotides)


    def set_format_key_value(self, key, value):
        '''Add a new key/value pair. Key in column 9 (FORMAT)
        and value in column 10. If key already exists, then updates
        the value to the new given value'''
        if self.format_keys is None:
            self.format_keys = []
            self.FORMAT = {}
        if key not in self.FORMAT:
            self.format_keys.append(key)
        self.FORMAT[key] = value


    def intersects(self, other):
        '''Returns True iff this record's reference positions overlap
        the other record reference positions (and are on same chromosome)'''
        return self.CHROM == other.CHROM and self.POS <= other.ref_end_pos() and other.POS <= self.ref_end_pos()


    def merge(self, other, reference_seq):
        '''Tries to merge this VcfRecord with other VcfRecord.
        Simple example (working in 0-based coords):
        ref = ACGT
        var1 = SNP at position 1, C->G
        var2 = SNP at position 3, T->A
        then this returns new variant, position=1, REF=CGT, ALT=GGA.

        If there is any kind of conflict, eg two SNPs in same position, then
        returns None.
        Also assumes there is only one ALT, otherwise returns None.'''
        if self.CHROM != other.CHROM or self.intersects(other) or len(self.ALT) != 1 or len(other.ALT) != 1:
            return None

        ref_start = min(self.POS, other.POS)
        ref_end = max(self.ref_end_pos(), other.ref_end_pos())
        ref_seq_for_vcf = reference_seq[ref_start:ref_end + 1]
        sorted_records = sorted([self, other], key=operator.attrgetter('POS'))
        alt_seq = []
        current_ref_pos = ref_start

        for record in sorted_records:
            assert record.REF != '.' and record.ALT[0] != '.'
            alt_seq.append(reference_seq[current_ref_pos:record.POS])
            alt_seq.append(record.ALT[0])
            current_ref_pos += len(record.REF)


        return VcfRecord('\t'.join([
            self.CHROM,
            str(ref_start + 1),
            '.',
            ref_seq_for_vcf,
            ''.join(alt_seq),
            '.', '.', 'SVTYPE=MERGED',
            'GT', '1/1',
        ]))

    def gt_aware_merge(self, other, reference_seq):
        '''Tries to merge this VcfRecord with other VcfRecord always using called allele as alt.
        Simple example (working in 0-based coords):
        ref = ACGT
        var1 = SNP at position 1, C->G called alt
        var2 = SNP at position 3, T->A called ref
        then this returns new variant, position=1, REF=CGT, ALT=GGT.

        If there is any kind of conflict, eg two SNPs in same position, then
        returns None.
        Also assumes there is only one ALT, otherwise returns None.'''
        if self.CHROM != other.CHROM or self.intersects(other) or len(self.ALT) != 1 or len(other.ALT) != 1:
            return None

        ref_start = min(self.POS, other.POS)
        ref_end = max(self.ref_end_pos(), other.ref_end_pos())
        ref_seq_for_vcf = reference_seq[ref_start:ref_end + 1]
        sorted_records = sorted([self, other], key=operator.attrgetter('POS'))
        alt_seq = []
        all_alt_seq = []
        gt_confs = []
        current_ref_pos = ref_start

        for record in sorted_records:
            assert record.REF != '.' and record.ALT[0] != '.'
            alt_seq.append(reference_seq[current_ref_pos:record.POS])
            all_alt_seq.append(reference_seq[current_ref_pos:record.POS])
            if record.FORMAT is None or 'GT' not in record.FORMAT:
                return None

            called_alleles = list(set(record.FORMAT['GT'].split('/')))
            if len(called_alleles) != 1 or '.' in called_alleles:
                return None
            gt = int(called_alleles[0])
            if gt > 0:
                alt_seq.append(record.ALT[gt-1])
            else:
                alt_seq.append(record.REF)
            all_alt_seq.append(record.ALT[0])
            current_ref_pos += len(record.REF)
            if record.FORMAT is not None and 'GT_CONF' in record.FORMAT:
                confs.append(record.FORMAT['GT_CONF'])

        alt_seq_for_vcf = ''.join(alt_seq)
        gt_conf = 0
        format = "GT"
        gt_0 = '0/0'
        gt_1 = '1/1'
        if len(gt_confs) > 0:
            gt_conf = min(gt_confs)
            format = 'GT:GT_CONF'
            gt_0 = '0/0:' + str(gt_conf)
            gt_1 = '1/1:' + str(gt_conf)

        if ref_seq_for_vcf == alt_seq_for_vcf:
            return VcfRecord('\t'.join([
                self.CHROM,
                str(ref_start + 1),
                '.',
                ref_seq_for_vcf,
                ''.join(all_alt_seq),
                '.', '.', 'SVTYPE=MERGED',
                format, gt_0,
            ]))
        else:
            return VcfRecord('\t'.join([
                self.CHROM,
                str(ref_start + 1),
                '.',
                ref_seq_for_vcf,
                alt_seq_for_vcf,
                '.', '.', 'SVTYPE=MERGED',
                format, gt_1,
            ]))


    def add_flanking_seqs(self, ref_seq, new_start, new_end):
        '''Adds new_start many nucleotides at the start, and new_end many nucleotides
        at the end from the appropriate nucleotides in reference sequence ref_seq.'''
        if new_start > self.POS or new_end < self.ref_end_pos():
            raise Error('new start and end positions must not try to shrink VCF record. new_start=' + str(new_start) + ', new_end=' + str(new_end) + '. VCF=' + str(self))

        new_start_nucleotides = ref_seq[new_start:self.POS]
        new_end_nucleotodes = ref_seq[self.ref_end_pos() + 1:new_end + 1]
        self.POS = new_start
        self.REF = new_start_nucleotides + self.REF + new_end_nucleotodes
        self.ALT = [new_start_nucleotides + x + new_end_nucleotodes for x in self.ALT]


    def merge_by_adding_new_alts(self, other, ref_seq):
        '''Adds other VcfRecord to this one, by adding new field(s) in the ALT
        column. Also adds REF nucleotides if they are needed. eg:
        ref: ACGT
        this var: pos=3, REF=T, ALT=A
        other var: pos=1, REF=C, ALT=T
        will change this var to be pos=1, REF=CGT, ALT=TGT,CGA'''
        if self.CHROM != other.CHROM:
            raise Error('Cannot merge two VCF records that lie on difference chromosomes\n' + str(self) + '\n' + str(other))

        new_ref_start = min(self.POS, other.POS)
        new_ref_end = max(self.ref_end_pos(), other.ref_end_pos())
        new_ref_end = min(new_ref_end, len(ref_seq) - 1)
        other = copy.copy(other)
        self.add_flanking_seqs(ref_seq, new_ref_start, new_ref_end)
        other.add_flanking_seqs(ref_seq, new_ref_start, new_ref_end)

        for other_alt in other.ALT:
            if other_alt not in self.ALT:
                self.ALT.append(other_alt)


    def ref_end_pos(self):
        '''Returns (zero-based) ref coord of the end of the variant'''
        return self.POS + len(self.REF) - 1


    def remove_useless_start_nucleotides(self):
        '''Removes duplicated nucleotides at the start of REF and ALT.
        But always leaves at least one nucleotide in each of REF and ALT.
        eg if variant is at position 42, REF=GCTGA, ALT=GCA, then
        sets position=41, REF=CTGA, ALT=CA.
        Assumes only one ALT, and does nothing if there is >1 ALT'''
        if len(self.REF) == 1 or len(self.ALT) != 1:
            return

        i = 0
        while i < len(self.REF) and i < len(self.ALT[0]) and self.REF[i] == self.ALT[0][i]:
            i += 1

        if i > 0:
            self.REF = self.REF[i - 1:]
            self.ALT = [self.ALT[0][i - 1:]]
            self.POS += i - 1


    def near_to_position(self, position, max_distance):
        '''Returns true iff the record is within max_distance of the given position.
        Note: chromosome name not checked, so that's up to you to do first.'''
        end = self.ref_end_pos()
        return self.POS <= position <= end or abs(position - self.POS) <= max_distance or abs(position - end) <= max_distance


    def inferred_var_seqs_plus_flanks(self, ref_seq, flank_length):
        '''Returns start position of first flank sequence, plus a list of sequences -
        the REF, plus one for each ALT.sequence. Order same as in ALT column'''
        flank_start = max(0, self.POS - flank_length)
        flank_end = min(len(ref_seq) - 1, self.ref_end_pos() + flank_length)
        seqs = [ref_seq[flank_start:self.POS] + self.REF + ref_seq[self.ref_end_pos() + 1: flank_end + 1]]

        for alt in self.ALT:
            seqs.append(ref_seq[flank_start:self.POS] + alt + ref_seq[self.ref_end_pos() + 1: flank_end + 1])

        return flank_start, seqs


    def total_coverage(self):
        '''Returns the sum of COV data, if present. Otherwise returns None'''
        if 'COV' in self.FORMAT:
            return sum([int(x) for x in self.FORMAT['COV'].split(',')])
        else:
            return None


    def called_alts_from_genotype(self):
        '''Returns a set of the (maybe REF and) ALT strings that were called, using GT in FORMAT.
        Returns None if GT not in the record'''
        if 'GT' not in self.FORMAT:
            return None

        genotype_indexes = set([int(x) for x in self.FORMAT['GT'].split('/')])
        alts = set()

        for i in genotype_indexes:
            if i == 0:
                alts.add(self.REF)
            else:
                alts.add(self.ALT[i-1])

        return alts


    def is_the_same_indel(self, other_record, ref_seq):
        '''Returns True iff this record and other_record are the "same"
        indel. At repeats, there is more than one way to report the same
        variant. eg:
        pos=42, ref=CAAA, alt=CAA
        pos=43, ref=AAA, alt=AA
        pos=44, ref=AA, alt=A'''
        if self.CHROM != other_record.CHROM or len(self.ALT) > 1 or len(other_record.ALT) > 1 or self.is_snp() or other_record.is_snp():
            return False

        # The number of nuleotides that have been added or removed
        # is a necessary condition of the indels being the same,
        # so check that before devling into the actual sequences
        if (len(self.REF) - len(self.ALT[0])) != (len(other_record.REF) - len(other_record.ALT[0])):
            return False

        # make records that start and end in the same place.
        # Then compare the REF and ALT sequences
        record1 = copy.copy(self)
        record2 = copy.copy(other_record)
        new_start = min(self.POS, other_record.POS)
        new_end = max(self.ref_end_pos(), other_record.ref_end_pos())
        record1.add_flanking_seqs(ref_seq, new_start, new_end)
        record2.add_flanking_seqs(ref_seq, new_start, new_end)
        return record1.REF == record2.REF and record1.ALT == record2.ALT


    def to_record_per_alt(self):
        '''Returns list of vcf_records. One per variant
        in the ALT column. Does not change INFO/FORMAT etc columns, which
        means that they are now broken'''
        record_list = []
        for alt in self.ALT:
            record_list.append(copy.copy(self))
            record_list[-1].ALT = [alt]
        return record_list


    def split_into_snps(self):
        '''Returns list of vcf_records. Tries to split
        this record into separate SNPs. eg if
        REF=ACGT and ALT=AGGA, then two SNPs
        C->G and T->A. Throws away all information in the
        INFO and FORMAT fields, except outputs the
        correct genotype (GT) if present in the input'''
        allele_lengths = set([len(x) for x in self.ALT])
        allele_lengths.add(len(self.REF))

        if len(allele_lengths) > 1 or allele_lengths == {1}:
            return [self]

        if self.FORMAT is not None and 'GT' in self.FORMAT:
            has_gt = True
            genotype_alleles = set([int(x) for x in self.FORMAT['GT'].split('/')])
        else:
            has_gt = False

        new_snps = {}

        for allele_index, allele in enumerate(self.ALT):
            for i in range(len(self.REF)):
                if self.REF[i] != allele[i]:
                    if i not in new_snps:
                        new_snps[i] = {'ref': self.REF[i], 'alts': {}}
                    assert new_snps[i]['ref'] == self.REF[i]
                    if allele[i] not in new_snps[i]['alts']:
                        new_snps[i]['alts'][allele[i]] = set()

                    if has_gt:
                        new_snps[i]['alts'][allele[i]].update(genotype_alleles.intersection({allele_index + 1}))


        new_vcfs = []

        for position_in_REF, allele_dict in sorted(new_snps.items()):
            new_vcfs.append(VcfRecord('\t'.join([
                self.CHROM,
                str(self.POS + position_in_REF + 1),
                '.',
                allele_dict['ref'],
                ','.join(sorted(list(allele_dict['alts'].keys()))),
                '.',
                'PASS',
                'SVTYPE=SNP',
            ])))

            if has_gt:
                if genotype_alleles == {0}:
                    gt = '0/0'
                else:
                    x = [len(allele_dict['alts'][x]) for x in sorted(allele_dict['alts'])]
                    matching_alleles = set([i+1 for i in range(len(x)) if x[i] > 0])

                    if len(matching_alleles) == 0:
                        gt = '0/0'
                    elif len(matching_alleles) == 1:
                        allele = matching_alleles.pop()
                        if len(genotype_alleles) == 1:
                            gt = str(allele) + '/' + str(allele)
                        else:
                            gt = '0/' + str(allele)
                    else:
                        assert len(matching_alleles) == 2
                        gt = '/'.join(sorted([str(x) for x in matching_alleles]))

                new_vcfs[-1].set_format_key_value('GT', gt)

        return new_vcfs

