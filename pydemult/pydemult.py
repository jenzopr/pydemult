#!/usr/bin/env python
import argparse
import logging
import itertools
import re
import multiprocessing
import subprocess
import pandas as pd
from mputil import lazy_map
from functools import partial

from .mutationhash import mutationhash
from .buffered_reader import buffered_blob
from .worker import _demult_chunk, _writer

# https://stackoverflow.com/a/43922107/6198494
def chunker_list(seq, size):
    return (seq[i::size] for i in range(size))

def demultiplex():
    parser = argparse.ArgumentParser(description='Demultiplexing of fastq files')
    parser.add_argument('--fastq', '-f', help='FASTQ file for demultiplexing.', metavar='input.fastq.gz', type=str)
    parser.add_argument('--samplesheet', '-w', help = 'Samplesheet containing barcodes and samplenames', metavar = 'samplesheet.txt', type=str)
    parser.add_argument('--barcode-regex', help = 'Regular expression to parse cell barcode (CB) and UMIs (UMI) from read names', metavar = '', default = 'none', type = str)
    parser.add_argument('--barcode-length', help='Length of cell barcode. Can be inserted into barcode regex via bclen.', metavar='11', type=str, default='11')
    parser.add_argument('--edit-distance', help='Maximum allowed edit distance for barcodes', metavar = '1', type=int, default = 1)
    parser.add_argument('--edit-alphabet', help='The alphabet that is used to created edited barcodes', choices=['N', 'ACGT', 'ACGTN'], default = "ACGT", type = str, metavar = "ACGT")
    parser.add_argument('--write-unmatched', help='Write reads with unmatched barcodes into unmatched.fastq.gz', action='store_true')
    parser.add_argument('--barcode-column', help='Name of the column containing barcodes', type=str, default='Barcode', metavar = 'Barcode')
    parser.add_argument('--sample-column', help='Name of the column containing sample names', type=str, default='Sample', metavar = 'Sample')
    parser.add_argument('--column-separator', help='Separator that is used in samplesheet', type=str, default='\t')
    parser.add_argument('--buffer-size', help="Buffer size for the FASTQ reader (in Bytes). Must be large enough to contain the largest entry.", type = int, default = 4000000, metavar = '4000000')
    parser.add_argument('--threads', '-t', help='Number of threads to use for multiprocessing.', type=int, metavar='1', default=1)
    parser.add_argument('--writer-threads', help='Number of threads to use for writing', type=int, metavar='1', default=2)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.3')
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
    # Create regular expression for barcode parsing from sequence header
    #
    if args.barcode_regex == "none":
        args.barcode_regex = '(.*):(?P<CB>[ATGCN]{#bclen#}'
    
    barcode_regex = args.barcode_regex.replace('#bclen#', args.barcode_length)
    #barcode_regex = '(.*):CELL_(?P<CB>[ATGCN]{'+str(args.barcode_length)+'}):UMI_(?P<UMI>[ATGCN]{8}):(.*)'
    #barcode_regex = '(.*):[ATGC]{0,5}(?P<CB1>[ATGCN]{6})TAGCCATCGCATTGC(?P<CB2>[ATGCN]{6})TACCTCTGAGCTGAA(?P<CB3>[ATGCN]{6})ACG(?P<UMI>[ATGCN]{6})GACT'
    c_barcode_regex = re.compile(barcode_regex)

    # Validate regex for presence of CB group and UMI group
    if 'CB' not in c_barcode_regex.groupindex.keys():
        logger.error('No cell barcode group CB found in barcode regex.')
        sys.exit(1)

    logger.debug('Set barcode regex to {}'.format(barcode_regex))

    #
    # Create mutationhash in case of barcode whitelist
    #
    samplesheet = pd.read_table(args.samplesheet, sep = args.column_separator, index_col = args.sample_column)
    barcode_dict = samplesheet[args.barcode_column].to_dict()
    barcodes = list(barcode_dict.values())

    logger.debug('Found barcodes:' + ','.join(barcodes))

    logger.debug('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.edit_distance, len(barcodes)))
    mut_hash = mutationhash(strings = barcodes, nedit = args.edit_distance, alphabet = list(args.edit_alphabet), log = logger)

    bufsize = args.buffer_size
    manager = multiprocessing.Manager()
    writer_pool = multiprocessing.Pool(args.writer_threads)

    # Handle writing queues
    if args.write_unmatched:
        queues = {'unmatched': manager.Queue()}
        queue_list = [queues['unmatched']]
        writer_pool.apply_async(_writer, (queues['unmatched'], ['unmatched']), callback = lambda x: print(x))
    else:
        queues = {}
        queue_list = []

    for chunk in chunker_list(list(barcode_dict.keys()), args.writer_threads-1):
        logger.debug('Creating writer queue for samples {}.'.format(','.join(chunk)))
        q = manager.Queue()
        q_bc_dict = dict((k, barcode_dict[k]) for k in chunk)
        writer_pool.apply_async(_writer, (q, q_bc_dict), callback = lambda x: print(x))
        queue_list.append(q)
        for bc in q_bc_dict.values():
            queues[bc] = q

    zcat = subprocess.Popen(['zcat', args.fastq], stdout=subprocess.PIPE, bufsize = bufsize)
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)
        demult_chunk = partial(_demult_chunk, mutationhash = mut_hash, regex = c_barcode_regex, write_unmatched = args.write_unmatched, q = queues)
        data = lazy_map(demult_chunk, blob_generator, n_cpus = args.threads)

    for q in queue_list:
        q.put((None, None))
    logger.debug('Shutting down queues.')

    writer_pool.close()
    writer_pool.join()
