import gzip
import re
import time
from functools import reduce

def entryfunc(blob):
    entries = blob.split(b'\n')
    # (header, sequence, quality)
    return(zip(entries[::4], entries[1::4], entries[3::4]))

def _demult_chunk(chunk, mutationhash, regex, write_unmatched, q):
    start = time.time()
    bc_chunks = dict()
    for entry in entryfunc(chunk):
        is_unmatched = False
        match = regex.match(entry[0].decode('utf-8'))

        if match is not None:
            try:
                bc_match = match.group('CB')
                origin = mutationhash[bc_match]

                if len(origin) > 1:
                    is_unmatched = True
                else:
                    barcode = list(origin)[0]
                    if barcode in bc_chunks:
                        bc_chunks[barcode].append(entry)
                    else:
                        bc_chunks[barcode] = [entry]
            except KeyError:
                is_unmatched = True

        if (match is None or is_unmatched) and write_unmatched:
            if 'unmatched' not in bc_chunks:
                bc_chunks['unmatched'] = [entry]
            else:
                bc_chunks['unmatched'].append(entry)
    parsing = time.time()
    for barcode in bc_chunks.keys():
        fastq = reduce(lambda x, y: x + '{}\n{}\n+\n{}\n'.format(y[0].decode('utf-8'), y[1].decode('utf-8'), y[2].decode('utf-8')), bc_chunks[barcode], '')
        q[barcode].put((barcode, fastq))
    writing = time.time()
    return((len(chunk), parsing-start, writing-parsing))

def _writer(q, barcodes):
    count = 0
    handles = dict()
    for (sample, barcode) in barcodes.items():
        handles[barcode] = gzip.open(sample + '.fastq.gz', 'wb')

    while 1:
        (bc, fastq) = q.get()
        if bc == None:
            break
        count += 1
        handles[bc].write(fastq.encode('utf-8'))

    for (sample, barcode) in barcodes.items():
        handles[barcode].close()
    return(count)
