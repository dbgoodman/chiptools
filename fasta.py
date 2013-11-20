import sys
import glob
import re
import os
import random
import warnings

from collections import defaultdict, namedtuple

from Bio import pairwise2, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_between_restr_sites(sr, rs1, rs2, prefix='', suffix='',
                                **kwargs):
    ''' 
    extract between two sites (i.e. restriction sites) and optionally append
    a prefix and suffix sequences (e.g. for adding illumina primers/barcodes)
    '''
    
    if isinstance(sr,list):
        return [extract_between_restr_sites(sr_i, rs1, rs2, prefix, suffix, 
                **kwargs) for sr_i in sr]
    elif isinstance(sr,SeqRecord):
        seq = sr.seq
    elif isinstance(sr,Seq):
        seq = sr
    else:
        raise(InputError,
            'Seq must be a SeqRecord list, a SeqRecord, or a Seq')
            
    if rs1 in seq.upper():
        index1 = seq.upper().find(rs1)
        if rs2 in seq.upper()[index1:]:
            index2 = seq.upper().find(rs2, index1)                
            
            extracted_seq = seq[index1+len(rs1):index2].tostring()
            appended_seq = ''.join([prefix, extracted_seq, suffix])
            
            if isinstance(sr,Seq):
                seq_id = ''
                seq_desc = ''
            elif isinstance(sr,SeqRecord):
                seq_id = sr.id
                seq_desc = sr.description
            return SeqRecord(seq=Seq(appended_seq),
                id=seq_id,
                description=seq_id)

        else: warnings.warn(
            'rs2 (%s) not found in sequence %s' % (rs2, seq))
    else: warnings.warn(
            'rs1 (%s) not found in sequence %s' % (rs1, seq))

def load_fasta_dir(path, **kwargs):
    '''
    load a directory full of .seq files, each containing one fasta record,
    into a seqrecord list
    '''
    seqrecords = []
    for infile in glob.glob(os.path.join(path, '*.seq') ):
      for record in SeqIO.parse(infile, 'fasta'):
          seqrecords.append(record)
    return seqrecords
 
def load_fasta_file(fasta_file, **kwargs):
   if isinstance(fasta_file,str):
       fasta_file = open(fasta_file)
      
   seqrecords = list(SeqIO.parse(fasta_file, "fasta"))
   return seqrecords

def gen_seq_unif(seqrecords):
    '''
    make a 'factory' object that randomly spits out fasta sequences
    '''
    seq_len = len(seqrecords)
    while True:
        i = random.randint(0,seq_len-1)
        yield (seqrecords[i].name, seqrecords[i].seq)

def clean_seq_name(sr, regex_str=r'^[A-Za-z-]*(\d+-*\d*)-(\d+).*', **kwargs):
    sr.description = sr.id
    (plate, num) = re.match(regex_str,sr.id).groups()
    sr.id = '%s-%d' % (plate, int(num))
    
def extract_sanger(seq_path, rs1, rs2,
                   output=None,
                   verbose=None,
                   revcom=False,
                   **kwargs):    
    print >> sys.stderr, 'Trimming and Extracting Sanger Seqs...'
    
    if os.path.isdir(seq_path):
        seqrecords = load_fasta_dir(seq_path, **kwargs)
    else:
        seqrecords = load_fasta_file(seq_path, **kwargs)
    
    [clean_seq_name(seq, **kwargs) for seq in seqrecords]

    if rs1 or rs2:
        seqrecords = extract_between_restr_sites(seqrecords, rs1, rs2)
    
    if revcom: 
        for sr in seqrecords: sr.seq = sr.seq.reverse_complement()
    
    if output:
        SeqIO.write(seqrecords,output,'fasta')
    else: return seqrecords
    
def extract_library(input_filename, rs1='GGCGCGCC', rs2='CATATG', 
                    output=None,
                    verbose=None,
                    **kwargs):
    print >> sys.stderr, 'Trimming and Extracting Library...'
    seqrecords = load_fasta_file(input_filename, **kwargs)        
    seqrecords = extract_between_restr_sites(seqrecords, rs1, rs2)
    if output:
        SeqIO.write(seqrecords,output,'fasta')
    else: 
        return seqrecords

