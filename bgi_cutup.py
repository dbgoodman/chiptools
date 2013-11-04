import os
import tempfile
import subprocess
import re
import chunk_fasta
import multiprocessing
import glob
import gzip
import sys
import time

raw_dir = '/scratch/tmp/bgihiseq/'

samples = {'ACCTGAAT':'202-1_R%d.fq.gz',
           'CAAGGCAT':'202-2_R%d.fq.gz', 
           'GTTCATAT':'203-1_R%d.fq.gz',
           'TATAAGAT':'203-2_R%d.fq.gz'}
           
filenames = []

def split_reads(fn, read):
    
    infile = gzip.open(fn)
    
    fobjs = {}
    filenames = []
    missing_total = 0
    read_total = 0
    
    #make files
    for fn in samples.values():
        tmp_file = os.path.join(raw_dir, fn % read)
        print "Making file "+tmp_file
        fobjs[fn % read] = gzip.open(tmp_file, 'w')
        filenames.append(tmp_file)
    
    try:
        while 1: 
            header = infile.next()
            read_total += 1
            bc = header[-11:-3]
            if bc not in samples.keys(): 
                missing_total += 1
                continue
                
            lines = "".join([header, 
                infile.next(), 
                infile.next(), 
                infile.next()])
            fobjs[samples[bc] % read].write(lines)
    except StopIteration:
        pass
    
    print "FILE "+infile.name+" COMPLETE: %d / %d: %f PCT BAD READS" % (missing_total, read_total, (missing_total/read_total*100))
    return filenames

def split_file(fn1, fn2, fileout, dir, read_num, bin):
    
    fileout.append('_%02d.fastq')
    fileout.append('.%02d.gz' % bin)
    print fn
    print ''.join(fileout)
    fileout = ''.join(fileout)
    
    lines = 0;
    try:
        while 1:
            header = infile.next()
            read_total += 1
            bc = header[-11:-3]
            if bc not in samples.keys(): 
                missing_total += 1
                continue
                
            lines = "".join([header, 
                infile.next(), 
                infile.next(), 
                infile.next()])
            print >> fobjs[samples[bc] % read], lines            
    gz.open
    
    
    chunk_fasta.store_chunk_files(fn, fileout, dir, 
         split_type= 'size', mb=250, mp= False)
    return True

# create a multiprocessing pool
procs = multiprocessing.cpu_count() - 1
pool = multiprocessing.Pool(procs)

jobs = list()

# jobs.append(pool.apply_async(split_reads, 
#     (raw_dir+"R1_concat.gz", 1), callback=filenames.extend))
# jobs.append(pool.apply_async(split_reads, 
#     (raw_dir+"R2_concat.gz", 2), callback=filenames.extend))
# while True:
#     done = sum([j.ready() for j in jobs])
#     tot = 2
# 
#     print >> sys.stderr, '\t%d / %d jobs completed.\r' % (
#              done, tot)
#     sys.stderr.flush()
# 
#     #if all tiles are complete:
#     if done == tot:
#         for job in jobs:
#             try: assert job.successful() 
#             except AssertionError: job.get()
#         break
#     else:
#         time.sleep(1)

filenames = ["/scratch/tmp/bgihiseq/202-1_R2.fq.gz",
    "/scratch/tmp/bgihiseq/203-2_R1.fq.gz",
    "/scratch/tmp/bgihiseq/202-2_R2.fq.gz",
    "/scratch/tmp/bgihiseq/202-2_R1.fq.gz",
    "/scratch/tmp/bgihiseq/203-1_R2.fq.gz",
    "/scratch/tmp/bgihiseq/203-1_R1.fq.gz",
    "/scratch/tmp/bgihiseq/203-2_R2.fq.gz",
    "/scratch/tmp/bgihiseq/202-1_R1.fq.gz"]

for fn in filenames:
        
    if 'R1' in fn:
        read_num = 1
    elif 'R2' in fn:
        read_num = 2
    
    if '-1' in fn:
        bin = 1
    elif '-2' in fn:
        bin = 2
    
    if '202' in fn:
        fileout = ['/scratch/dbg/ecre/fq/202_hsrna/s_G1_L001_%s' % read_num]
        dir = '/scratch/dbg/ecre/fq/202_hsrna/'
    elif '203' in fn:
        fileout = ['/scratch/dbg/ecre/fq/203_hsrna/s_G1_L001_%s' % read_num]
        dir = '/scratch/dbg/ecre/fq/203_hsrna/'
        
    print fileout
    
    pool.apply_async(split_file, [fn, fileout, dir, read_num, bin])
    #split_file(fn, fileout, dir, read_num, bin)

# close up the pool if we no longer want to swim
pool.close()
pool.join()
        