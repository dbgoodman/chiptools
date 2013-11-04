import itertools
import pandas as p
import numpy as np
import subprocess

from collections import defaultdict

#total RNA reads in bin (tR)
#total DNA reads in bin (tD)

# RNA:
# for each library member:
#     - how many RNA members match
#     - total RNA read count
#     - total DNA read count
#     - how many RNA members match uniquely
#     - unique RNA read count
#     - total RNA / tR : DNA ratio / tD
#     - unique RNA : DNA ratio
#     - percentage of RNA members that map uniquely / total RNA members
#     - percentage of read counts that map uniquely / total read count
#     - TSS positions per unique RNA
#     - total read counts per TSS position

def load_fastq(fastq_fn):
    
    fq_seqs = dict()
    
    with open(fastq_fn) as f:
        while f:
            try:
                seq_id = f.next().lstrip('>').rstrip()
                seq_name = f.next().rstrip()
                fq_seqs[seq_id] = seq_name
            except(StopIteration):
                f.close()
                break        
    
    return fq_seqs

def run_bowtie(lib_prefix, prefix_fn, 
    unique_fa, bowtie_cmd= None, min_read_count = 5, max_read_length = 200):
    
    if not bowtie_cmd:
        # remove the newlines from all headers and change them to tabs 
        # then feed the two-column file to bowtie, then filter on minimum
        # read count and maximum read length through perl
        bowtie_cmd = '''
            bowtie -v 3 -l 10 -k 114 -p 16 \
                --norc --best --strata --suppress 2,6 \
                --un %(prefix_fn)s.unmapped.fa -f %(lib_prefix)s  \
                <(perl -pe 's/([^ATGC])\\n/$1\\t/' %(unique_fa)s \
                | perl -ne '@l = split; ($l[1] > %(min_read_count)d 
                    && length($l[4]) < %(max_read_length)d) 
                    && (s/\\t([ATGC])/\\n$1/ && print);')''' % locals()
                        
    p = subprocess.Popen(bowtie_cmd % locals(), 
        stdout= subprocess.PIPE,
        shell= True,
        executable='/bin/bash')
        
    return p.stdout

def total_fasta_reads(fasta_fn):
    
    total = None
    
    with open(fasta_fn) as f:
        while f:
            try:
                l = f.next()
                if l.startswith('>'):
                    if not total: 
                        total = [0] * len(l.strip().split()[2:])
                    total = [int(n) + total[i] for i,n in 
                        enumerate(l.strip().split()[2:])]
            except(StopIteration):
                f.close()
                break
    return total

def load_aligns(lib_fasta_fn, align_fasta_fn, type):
    ''' 
    This file will be in bowtie format:
     http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output
     
    My read names are 4 tab-separated fields, read number, count, bin 1, bin 2.
    
    01 Name of read that aligned
    02 Total read count
    03 Read count in bin 1
    04 Read count in bin 2
    (suppressed) Reference strand aligned to
    05 Name of reference sequence where alignment occurs
    06 0-based offset into the forward reference strand 
    07 Read sequence
    (suppressed) Read qualities
    08 Number of other alignments for this RNA
    09 Mismatches (base:N>N, ... )
    '''
    
    #load library fasta file into a dict
    lib_seqs = load_fastq(lib_fasta_fn)
    
    #bowtie fields
    field_names = ['read_num', 'tot_count', 'rb1_count', 'rb2_count',
        'lib_id', 'l_offset', 'seq', 'alt_aligns', 'mismatches']
    field_dtypes = "i8, i8, i8, i8, S40, i8, S200, i8, S40"
    
    #run bowtie to align reads
    bowtie_output = run_bowtie(lib_fasta_fn.rstrip('.fa'), prefix_fn+type, 
        align_fasta_fn)
    rna_rec = p.DataFrame(np.genfromtxt(
        fname= bowtie_output, 
        dtype=field_dtypes,
        delimiter="\t",
        names=field_names))


lib_fasta_fn = "/scratch/dbg/ecre/fa/203.norestrict.fa"
rna_fasta_fn = "/scratch/dbg/ecre/ct/203_hsrna/203_hsrna.counts.fa"
dna_fasta_fn = "/scratch/dbg/ecre/ct/203_hsdna/203_hsdna.counts.fa"
prefix_fn = "/scratch/dbg/ecre/203_hs"

# load both RNA and DNA
rna_rec = load_aligns(lib_fasta_fn, rna_fasta_fn, type='rna')
dna_rec = load_aligns(lib_fasta_fn, dna_fasta_fn, type='dna')

# RNA PROCESSING===============================================================
# How do I determine if I want to keep an RNA?
# 1. does it map without any mismatches?
# 2. does it map all the way to the end of a library member?
# 3. does it map uniquely to this library member?
# 4. does it map in the middle (i.e. not DNA)



rna_rec = p.concat([
    rna_rec, 
    p.DataFrame([dict(zip( ('promoter','rbs'), i.split('--'))) 
        for i in  rna_rec['lib_id']]
    )], axis=1)

rna_rec_uniq = (rna_rec[
    (rna_rec['r_offset'] == 0) 
    #& ( (rna_rec['tot_count'] > 100) | (rna_rec['mismatches'] == '') )
    #the BBa_J61100 hack is b/c Anderson matches with r_offset == 1 \/
    & ( (rna_rec['alt_aligns'] == 0) | 
        ((rna_rec['rbs'] == 'BBa_J61100') & (rna_rec['alt_aligns'] == 1)) )
    & (rna_rec['l_offset'] > 2)])

rna_rec_not_used = (rna_rec[
    ~( (rna_rec['r_offset'] == 0) 
    #& ( (rna_rec['tot_count'] > 100) | (rna_rec['mismatches'] == '') )
    #the BBa_J61100 hack is b/c Anderson matches with r_offset == 1 \/
    & ( (rna_rec['alt_aligns'] == 0) | 
        ((rna_rec['rbs'] == 'BBa_J61100') & (rna_rec['alt_aligns'] == 1)) )
    & (rna_rec['l_offset'] > 2) ) & (rna_rec['tot_count'] >= 500)])

#print out missing 203 RNA members
rna_lib_missing = (set(load_fastq(lib_fasta_fn).keys()) - 
    set(p.unique(rna_rec_uniq['lib_id'])))
    
# missing203rna = set([item for item in rna_lib_missing if not '--' in item])
# print "%d RNA sequences from 203 are missing." % len(missing203rna)

#save a list of raw read data
rna_rec.to_csv(prefix_fn+'.raw_rna.csv')
rna_rec_uniq.to_csv(prefix_fn+'.uniq_rna.csv')

rna_rec_not_used.to_csv(prefix_fn+'.not_used.csv')

#save a summed list of read counts per bin
rna_rb_counts = rna_rec_uniq[
    ['lib_id','rb1_count','rb2_count']].groupby('lib_id').sum()

rna_rb_counts.to_csv(prefix_fn+'.rna_counts.csv')

best_record_per_group = lambda df: df.sort_index(
    ascending=False, by='tot_count').head(2)

best_records = rna_rec[
    (rna_rec['r_offset'] == 0) 
    & ( (rna_rec['tot_count'] > 100) | (rna_rec['mismatches'] == '') )
    & (rna_rec['alt_aligns'] < 40) & (rna_rec['l_offset'] > 2)
    ].groupby('lib_id').apply(best_record_per_group)

odd_records = best_records[(best_records['alt_aligns'] > 0)
    | (best_records['mismatches'] != '')]

# DNA PROCESSING==============================================================

dna_rec = p.concat([
    dna_rec, 
    p.DataFrame([dict(zip( ('promoter','rbs'), i.split('--'))) 
        for i in  dna_rec['lib_id']]
    )], axis=1)

dna_rec_uniq = (dna_rec[
    (dna_rec['r_offset'] == 0) 
    #& (dna_rec['mismatches'] == '') 
    & ((dna_rec['alt_aligns'] == 0) |
        ((dna_rec['rbs'] == 'BBa_J61100') & (dna_rec['alt_aligns'] == 1)) )
    & (dna_rec['l_offset'] == 0)
    & (dna_rec['rb1_count'] > 5)
    & (dna_rec['rb2_count'] > 5)])

#print out missing 203 DNA members
dna_lib_missing = (set(load_fastq(lib_fasta_fn).keys()) - 
    set(p.unique(dna_rec_uniq['lib_id'])))
    
# missing203dna = set([item for item in dna_lib_missing if not '--' in item])
# print "%d DNA sequences from 203 are missing." % len(missing203dna)

#save a summed list of read counts per bin
dna_rb_counts = dna_rec_uniq[
    ['lib_id','rb1_count','rb2_count']].groupby('lib_id').sum()

dna_rb_counts.to_csv(prefix_fn+'.dna_counts.csv')

# DNA/RNA COMPARISONS=========================================================

# rna_not_dna = missing203rna - missing203dna
# dna_not_rna = missing203dna - missing203rna

rna_not_dna = rna_lib_missing - dna_lib_missing
dna_not_rna = dna_lib_missing - rna_lib_missing


print "%d missing RNA sequences are present in DNA." % len(rna_not_dna)
print "%d missing DNA sequences are present in RNA." % len(dna_not_rna)

#HARD CODE THIS INSTEAD \/ (FROM UNIQ READ COUNTS BEFORE BOWTIE)
# rna_sums = rna_rec_uniq[['rb1_count','rb2_count']].sum()
# dna_sums = dna_rec_uniq[['rb1_count','rb2_count']].sum()
rna_sums = total_fasta_reads(rna_fasta_fn)
dna_sums = total_fasta_reads(dna_fasta_fn)

rna_ratio = (rna_rb_counts / rna_sums)
rna_avgs = (rna_ratio['rb1_count'] + rna_ratio['rb2_count']) / 2 

dna_ratio = (dna_rb_counts / dna_sums)
dna_avgs = (dna_ratio['rb1_count'] + dna_ratio['rb2_count']) / 2 

rna_to_dna_ratio = rna_avgs / dna_avgs

rna_to_dna_ratio.to_csv(prefix_fn+'.rna_dna_ratio.csv')
