#!/usr/bin/env python
import gzip
import bz2

def extract_nr(fields):
    return fields.split(':')[4]

def extract_nv(fields):
    return fields.split(':')[5]

def extract_sc(fields):
    for data in fields.split(';'):
        key, val = data.split('=')
        if key == 'SC':
            return val

def detect_compression(filename):
    with open(filename, 'rb') as fl:
        magic = fl.read(3)

    if magic[:2] == b'\x1f\x8b':
        return 'GZIP'

    elif magic[:3] == b'BZh':
        return 'BZIP2'

    else:
        return None

def run(in_filename, out_filename):
    def _do_work(ifl, ofl):
        for linenumber, line in enumerate(ifl, start=1):
                if linenumber % 10000 == 0:
                    print('processed {} lines'.format(linenumber))
                if line.startswith('##'):
                    continue

                elif line.startswith('#'):
                    colheaders = line.rstrip().split('\t')
                    sample_names = colheaders[9:]
                    nsamples = len(sample_names)
                    sample_nrs = ['{}.nr'.format(name) for name in sample_names]
                    sample_nvs = ['{}.nv'.format(name) for name in sample_names]
                    out_colheaders = [colheaders[0]]
                    for n in [1,3,4,5,7]:
                        out_colheaders.append(colheaders[n])
                    for n in sample_nrs + sample_nvs:
                        out_colheaders.append(n)
                    ofl.write('\t'.join(out_colheaders) + '\n')

                else:
                    if nsamples is None:
                        raise RuntimeError('Sample names must not be None - header was not parsed correctly')
                    fields = line.rstrip().split('\t')
                    out_fields = [fields[0]] # chrom
                    out_fields.append(fields[1]) # pos
                    out_fields.append(fields[3]) # from
                    out_fields.append(fields[4]) # to
                    out_fields.append(fields[5]) # qual
                    out_fields.append(extract_sc(fields[7])) # context
                    nrs = [extract_nr(fields[9+i]) for i in range(nsamples)]
                    nvs = [extract_nv(fields[9+i]) for i in range(nsamples)]
                    out_fields.extend(nrs)
                    out_fields.extend(nvs)
                    ofl.write('\t'.join(out_fields) + '\n')

    nsamples = None

    compression_type = detect_compression(in_filename)

    if compression_type == 'GZIP':
        with gzip.open(in_filename, 'rb') as ifl:
            with gzip.open(out_filename, 'wb') as ofl:
                _do_work(ifl, ofl)

    elif compression_type == 'BZIP2':
        with bz2.BZ2File(in_filename, 'rb') as ifl:
            with bz2.BZ2File(out_filename, 'wb') as ofl:
                _do_work(ifl, ofl)

    else:
        with open(in_filename, 'r') as ifl:
            with open(out_filename, 'w') as ofl:
                _do_work(ifl, ofl)


if __name__ == '__main__':
    import os, sys
    if len(sys.argv) < 3:
        print('process_vcf.py extracts read coverage and variant read counts from a somatypus VCF.')
        print('Also extracts sequence context, ref and alt base, genomic position and quality score')
        print('Input can be uncompressed, gzipped or bzip2ed. Output file compression will match input')
        print('Usage: (python) process_vcf.py INFILE OUTFILE')
        sys.exit(1)
    infile, outfile= sys.argv[1:3]
    if not os.path.isfile(infile):
        print('Input file {} not found'.format(infile))
        sys.exit(1)
    if not os.access(os.path.dirname(os.path.abspath(outfile)), os.W_OK):
        print('Output file {} is not in a writeable location'.format(outfile))
        sys.exit(1)
    if detect_compression(infile) == 'GZIP' and not outfile.endswith('.gz'):
        outfile = outfile+('.gz')
    elif detect_compression(infile) == 'BZIP2' and not outfile.endswith('.bz2'):
        outfile = outfile + ('.bz2')
    run(infile, outfile)
