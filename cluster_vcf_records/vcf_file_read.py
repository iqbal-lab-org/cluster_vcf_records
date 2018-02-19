import operator
import logging

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

