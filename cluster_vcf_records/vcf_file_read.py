import operator
import logging

import pyfastaq

from cluster_vcf_records import vcf_record

class Error (Exception): pass

def vcf_file_to_dict(infile, sort=True, homozygous_only=False, remove_asterisk_alts=False, max_REF_len=None, remove_useless_start_nucleotides=False):
    header_lines = []
    records = {}
    f = pyfastaq.utils.open_file_read(infile)
    count_keys = ['keep', 'homozygous', 'REF_too_long', 'alt_is_asterisk']
    counts = {x: 0 for x in count_keys}

    for line in f:
        if line.startswith('#'):
            header_lines.append(line.rstrip())
            continue

        record = vcf_record.VcfRecord(line)
        if homozygous_only and record.FORMAT.get('GT', None) != '1/1':
            counts['homozygous'] += 1
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

