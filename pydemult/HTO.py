import sys
import argparse
import logging
import re

from .mutationhash import mutationhash

def count():
    parser = argparse.ArgumentParser(description='Demultiplex samples based fastq files from hash tag oligo data')
    parser.add_argument('--reference', '-r', help='Tab-separated reference file containing hash tag sequences and names')
    parser.add_argument('--whitelist', '-w', help='Cell barcode whitelist of allowed / known cell barcodes')
    parser.add_argument('--barcode-regex', '-b', help = 'Regular expression to parse cell barcode (CB) from barcode sequences', default = '(.*):(?P<CB>[ATGCN]{11}', type = str)
    parser.add_argument('--hashtag-regex', '-h', help = 'Regular expression to parse hash tag sequences (HTO) from hash tag sequences', default = '(.*)(?P<HTO>[ATGCN]{15}', type = str)
    parser.add_argument('--barcode-edit-distance', help='Maximum allowed edit distance for barcodes', metavar = '1', type=int, default = 1)
    parser.add_argument('--hashtag-edit-distance', help='Maximum allowed edit distance for hash tag oligos', metavar = '2', type=int, default = 1)
    parser.add_argument('--edit-alphabet', help='The alphabet that is used to created edited barcodes / hash tag sequences', choices=['N', 'ACGT', 'ACGTN'], default = "ACGTN", type = str, metavar = "ACGTN")
    parser.add_argument('--buffer-size', help="Buffer size for the FASTQ reader (in Bytes). Must be large enough to contain the largest entry.", type = int, default = 4000000, metavar = '4000000')
    parser.add_argument('--threads', '-t', help='Number of threads to use for multiprocessing.', type=int, metavar='1', default=1)
    parser.add_argument('--barcode-sequences', '-i', help='FASTQ file containing barcode sequences', metavar='input_BC.fastq.gz', type=str)
    parser.add_argument('--hashtag-sequences', '-j', help='FASTQ file containing hash tag sequences', metavar='input_HT.fastq.gz', type=str)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(logging.DEBUG)
    logger.addHandler(sh)

    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)

    logger.debug('Working on {} using {} threads'.format(args.fastq, args.threads))

    #
    # Create regular expression for barcode / hash tag oligo parsing from sequences
    #
    c_barcode_regex = re.compile(args.barcode_regex)
    c_hashtag_regex = re.compile(args.hashtag_regex)

    # Validate regex for presence of CB group
    if 'CB' not in c_barcode_regex.groupindex.keys():
        logger.error('No cell barcode group CB found in barcode regex.')
        sys.exit(1)
    # Validate regex for presence of HTO group
    if 'CB' not in c_hashtag_regex.groupindex.keys():
        logger.error('No cell barcode group HTO found in hashtag regex.')
        sys.exit(1)

    #
    # Create mutationhash for barcodes from whitelist
    #
    with open(args.whitelist, 'r') as file:
        barcodes = [line.rstrip('\n') for line in file]

    logger.debug('Whitelist contains the following barcodes' + ','.join(barcodes))

    logger.debug('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.edit_distance, len(barcodes)))
    # mut_hash = mutationhash(strings = barcodes, nedit = args.edit_distance, alphabet = list(args.edit_alphabet), log = logger)

    
