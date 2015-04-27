import exceptions
import subprocess

import config

def num(s):
    '''convert a value to either a int or a float'''
    try:
        return int(s)
    except exceptions.ValueError:
        return float(s)

def open_zcat(fn):
        if not isinstance(fn,list):
            fn = [fn]

        p = subprocess.Popen([config.bin_paths['zcat_path']]+fn,
            stdout = subprocess.PIPE)

        return p.stdout

def quick_fasta(fa_file, get_annotations= False, seqdict = {}):
    '''quickly read from a fasta file without making objects in biopython'''

    if isinstance(fa_file, str):
        fa_fhandle = open(fa_file, 'r')
    else:
        fa_fhandle = fa_file

    proceed = True
    try:
        next_header = fa_fhandle.next().rstrip().split()
    except StopIteration:
        proceed = False

    while proceed:
        try:
            header = next_header
            seq = fa_fhandle.next().rstrip()
            next_header = fa_fhandle.next().rstrip().split()

            while not next_header[0].startswith('>'):
                seq += next_header[0]
                next_header = fa_fhandle.next().rstrip().split()

        except StopIteration:
            proceed = False
        finally:
            seqid = header[0][1:]

            if get_annotations:
                seqdict[seqid] = ( seq.rstrip(),
                    eval(' '.join(header[2:])) )
            else:
                seqdict[seqid] = (seq.rstrip(), )

    return seqdict
