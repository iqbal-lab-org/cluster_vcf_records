import itertools
import operator
import logging
import multiprocessing

import pyfastaq

from cluster_vcf_records import vcf_record

class Error (Exception): pass

def vcf_file_to_dict(infile, sort=True, homozygous_only=False, remove_asterisk_alts=False, max_REF_len=None, remove_useless_start_nucleotides=False, min_SNP_qual=None, min_dp4=None, min_GT_conf=None):
    header_lines = []
    records = {}
    f = pyfastaq.utils.open_file_read(infile)
    count_keys = ['keep', 'not_homozygous', 'REF_too_long', 'alt_is_asterisk','SNP_qual_too_low','not_enough_high_qual_reads','GT_conf_too_low','samtools_indel','cortex_snp_call']
    counts = {x: 0 for x in count_keys}

    for line in f:
        if line.startswith('#'):
            header_lines.append(line.rstrip())
            continue

        record = vcf_record.VcfRecord(line)
        if homozygous_only and record.FORMAT.get('GT', None) != '1/1':
            counts['not_homozygous'] += 1
            continue

        if min_GT_conf != None:
            if record.QUAL == None:

                if record.INFO.get('SVTYPE') == 'SNP':      ##Clears out any cortex SNP calls
                    counts['cortex_snp_call'] += 1
                    continue

                if float(record.FORMAT.get('GT_CONF')) < min_GT_conf:   ##Filters based on cortex gt_conf
                    counts['GT_conf_too_low'] += 1
                    continue

            if record.QUAL != None:

                #if 'INDEL' in record.INFO:        ##Clears out any samtools indel calls
                    #counts['samtools_indel'] += 1
                    #continue

                if min_SNP_qual != None and record.QUAL < min_SNP_qual:    ##Filters based on samtools qual score
                    counts['SNP_qual_too_low'] += 1
                    continue

                if record.INFO.get('DP4') != None:          ##Filters for at least [] dp4 alt reads on both strands
                    if min_dp4 != None:
                        dp4_cov = record.INFO['DP4']
                        cov_split = dp4_cov.split(',')
                        if float(cov_split[2]) < min_dp4 or float(cov_split[3]) < min_dp4:
                            counts['not_enough_high_qual_reads'] += 1
                            continue

        if remove_asterisk_alts:
            record.remove_asterisk_alts()

        if len(record.ALT) < 1:
            counts['alt_is_asterisk'] += 1
            continue

        if remove_useless_start_nucleotides:
            record.remove_useless_start_nucleotides()

        if max_REF_len is not None and len(record.REF) > max_REF_len:
            counts['REF_too_long'] += 1
            continue

        if record.CHROM not in records:
            records[record.CHROM] = []
        records[record.CHROM].append(record)
        counts['keep'] += 1

    pyfastaq.utils.close(f)
    logging.info('Loaded file ' + infile + '. Counts: ' + ';'.join([x + '=' + str(counts[x]) for x in count_keys]))

    if sort:
        for record_list in records.values():
            record_list.sort(key=operator.attrgetter('POS'))

    return header_lines, records


def vcf_file_to_list(infile):
    header_lines = []
    records = []
    f = pyfastaq.utils.open_file_read(infile)

    for line in f:
        if line.startswith('#'):
            header_lines.append(line.rstrip())
        else:
            records.append(vcf_record.VcfRecord(line))

    pyfastaq.utils.close(f)
    return header_lines, records


def get_sample_name_from_vcf_header_lines(header_lines):
    '''Given a list of header lines (made by either
    vcf_file_to_dict() or vcf_file_to_list()),
    returns the sample name. Assumes only one sample
    in the file.
    Raises error if badly formatted #CHROM line.
    Returns None if no #CHROM line found'''
    # We want the #CHROM line, which should be the last line
    # of the header
    for line in reversed(header_lines):
        if line.startswith('#CHROM'):
            fields = line.rstrip().split('\t')
            required_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            if fields[:len(required_cols)] != required_cols:
                raise Error('Error! #CHROM line must have these for first 8 columns: ' + ', '.join(required_cols) + '\nat this line of file: ' + line)

            if len(fields) == len(required_cols):
                return None

            required_cols.append('FORMAT')
            format_column = fields[len(required_cols) - 1]
            if format_column != required_cols[-1]:
                raise Error('Error! #CHROM line has 9^th column, which should be "FORMAT" but is: "' + format_column + '" at this line of file: ' + line)

            if len(fields) == len(required_cols):
                logging.warning('FORMAT column in file, but no sample names at line: ' + line.rstrip())
                return None

            sample = fields[len(required_cols)]

            if len(fields) > len(required_cols) + 1:
                logging.warning('More than one sample found. Using the name of the first sample (' + sample + ') at line: ' + line.rstrip())

            return sample

    return None


def vcf_file_has_at_least_one_record(infile):
    f = pyfastaq.utils.open_file_read(infile)

    for line in f:
        if line.startswith('#'):
            pass
        else:
            record = vcf_record.VcfRecord(line)
            pyfastaq.utils.close(f)
            return True

    pyfastaq.utils.close(f)
    return False


def get_header_lines_from_vcf_file(infile):
    f = pyfastaq.utils.open_file_read(infile)
    header_lines = [line.rstrip() for line in f if line.startswith('#')]
    pyfastaq.utils.close(f)
    return header_lines



def get_sample_name_from_vcf_file(infile):
    '''Assumes only one sample in the file.
    Raises error if badly formatted #CHROM line.
    Returns None if no #CHROM line found'''
    header_lines =  get_header_lines_from_vcf_file(infile)
    return get_sample_name_from_vcf_header_lines(header_lines)



def vcf_file_to_dict_of_vars(infile, reference_seqs):
    '''Loads just the variant info from input VCF file.
    If reference_seqs is given, should be a dict of seq name -> pyfastaq Fastaq sequence,
    and will be used to sanity check variants in input file. Any where CHROM
    is not in the dict, or REF string does not match the ref sequence, are ignored
    Output is a dictionary of:
    ref name -> position -> {ref string -> {set of alt strings}}'''
    variants = {}
    _, vcf_records = vcf_file_to_dict(infile, sort=False, remove_asterisk_alts=True, remove_useless_start_nucleotides=True)

    for ref_name, variant_list in vcf_records.items():
        for record in variant_list:
            if record.POS < 0:
                logging.warning(f'VCF record with negative POS in file {infile}. Ignoring: {record}')
                continue
            elif record.CHROM not in reference_seqs:
                logging.warning(f'CHROM not recognised in VCF record in file {infile}. Ignoring: {record}')
                continue
            elif reference_seqs[record.CHROM][record.POS:record.POS + len(record.REF)] != record.REF:
                logging.warning(f'REF string does not match reference seq in file {infile}. Ignoring: {record}')
                continue

            if ref_name not in variants:
                variants[ref_name] = {}

            if record.POS not in variants[ref_name]:
                variants[ref_name][record.POS] = {}

            if record.REF not in variants[ref_name][record.POS]:
                variants[ref_name][record.POS][record.REF] = set()

            variants[ref_name][record.POS][record.REF].update(record.ALT)

    return variants


def vcf_files_to_dict_of_vars(infiles, reference_seqs, threads=1):
    variants = {}

    for i in range(0, len(infiles), threads):
        with multiprocessing.Pool(threads) as pool:
            new_variants_dict_list = pool.starmap(vcf_file_to_dict_of_vars, zip(infiles[i:i+threads], itertools.repeat(reference_seqs)))

        for new_variants in new_variants_dict_list:
            for ref_name in new_variants:
                if ref_name not in variants:
                    variants[ref_name] = new_variants[ref_name]
                    continue

                for pos in new_variants[ref_name]:
                    if pos not in variants[ref_name]:
                        variants[ref_name][pos] = new_variants[ref_name][pos]
                        continue

                    for ref, new_alts in new_variants[ref_name][pos].items():
                        if ref not in variants[ref_name][pos]:
                            variants[ref_name][pos][ref] = new_alts
                        else:
                            variants[ref_name][pos][ref].update(new_alts)

        if i % 100 == 0:
            logging.info(f'Loaded {i+threads} files out of {len(infiles)}')

    return variants

