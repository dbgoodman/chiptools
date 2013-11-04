#!/usr/bin/env python
# encoding: utf-8

'''
Created by Brant Faircloth on 11 December 2010 11:28 PST (-0800).
MODIFIED HEAVILY BY ME, ADDED OPTIONS, SPLITTING, AND GZIP COMPATIBILITY
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

*   Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

*   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

*   Neither the name of the University of California nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission. THIS SOFTWARE IS
    PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    POSSIBILITY OF SUCH DAMAGE. 

'''

import os
import tempfile
import multiprocessing
import gzip
import shutil
import itertools

import pdb

def file_type(input):
    """given an input file, determine the type and return both type and record delimiter (> or @)"""
    name, extension = os.path.splitext(os.path.basename(input))
    fastas = set(['.fsa','.fasta','.fa'])
    fastqs = set(['.fastq','.fq'])
    fq_gz = set(['.gz','.gz'])
    gffs = set(['.gff'])
    gz = False
    if extension in fastas:
        ft = 'fasta'
        delim = '>'
    elif extension in fastqs:
        ft = 'fastq'
        delim = '@HWI-' #I added the HWI b/c some quality lines start with @!
    elif extension in fq_gz:
        ft = 'fastq'
        delim = '@HWI-' #I added the HWI b/c some quality lines start with @!
        gz = True
    else:
        raise IOError, "Input file not of correct extension"
    return ft, delim, gz
    
def _get_file_chunks(input, delim, size, gz):
    """given input, record delimiter, and chunk size, yield an iterator contains file seek (start)
    and file read (stop) positions.  Return final position as (6365605, None)."""
    if gz:
        f = gzip.open(input)
    else:
        f = open(input)
    while 1:
        start = f.tell()
        f.seek(size, 1)
        line = f.readline()
        if not line:
            break
        # if this isn't a fasta header line, read forward until
        # we get to one
        while not line.startswith(delim):
            line = f.readline()
        else:
            # now that we got to a fasta header, we're at the end.
            # back up the length of the fasta header.
            f.seek(-len(line), 1)
            # tuple up
            yield start, f.tell() - start, input
    # make sure we catch the (start, distance) for the end of the file, too
    yield start, None, input
    f.close()
    
def get_chunks(input, delim, split_type, mb=1, splits= None, gz= False):
    
    #do not split below 128 kb
    minsize = (1024**2) / 8
    
    """return a tuple of file seek (start, distance) positions covering chunks of a file"""
    if os.path.getsize(input) > minsize:
        if split_type == 'size':
            size = mb * (1024**2)
        if split_type == 'pieces':
            if not splits:
                splits = multiprocessing.cpu_count() - 1
            size = int(round((os.path.getsize(input)/float(splits)), 0))
    else:
        size = os.path.getsize(input)
    return _get_file_chunks(input, delim, size, gz)

def _split_file(chunk_gz):
    """function to split a file into appropriate pieces given (start, stop) file seek coords"""
    
    chunk, gz = chunk_gz
    
    if gz:
        f = gzip.open(chunk[2])
    else:
        f = open(chunk[2])
    f.seek(chunk[0])
    if chunk[1]:
        d = f.read(chunk[1])
    else:
        d = f.read()
    td, tf = tempfile.mkstemp(suffix='.splt')
    print "Making tempfile %s" % tf
    os.close(td)
    if gz:
        otf = gzip.open(tf, 'wb')
    else:
        otf = open(tf, 'w')
    print "\t" + tf
    otf.write(d)
    otf.close()
    f.close()
    return tf
    
def make_chunks(chunks, pool = None, mp = True, gz= False):
    """return a list of tempfiles that are the chunked input file"""
    chunk_gz_iter = itertools.izip(chunks, itertools.repeat(gz))
    
    mp = False
    
    if mp and not pool:
        # create a multiprocessing pool
        procs = multiprocessing.cpu_count() - 1
        pool = multiprocessing.Pool(procs)
        chunks = pool.map(_split_file, chunk_gz_iter)
        # close up the pool if we no longer want to swim
        pool.close()
        pool.join()
    else:
        chunks = map(_split_file, chunk_gz_iter)
        
    return chunks

def chunk_file(input_fn, split_type='pieces', mb= 24, splits= None, 
    tempdir=tempfile.tempdir, mp= True, **kwargs):
    
    #this needs to be set manually in some cases so that we can easily move 
    #file names after making temporary files without shuttling all the data 
    #between mounts
    tempfile.tempdir = tempdir
    
    f_type, delim, gz = file_type(input_fn)
    chunk_offsets = get_chunks(input_fn, delim, split_type, mb, splits, gz)
    chunks = make_chunks(chunk_offsets, mp= mp, gz= gz)
    return chunks

def store_chunk_files(input_fn, output_form, output_dir, **kwargs):
    '''In case you want your chunks to be a little less temporary, but also
    to be gzipped, this allows you to do so given an output format. It
    makes temporary chunk files and then changes their filenames according 
    to the format specified. 
    
    output form should have one format digit for the chunk number. For
    instance, if the input file is:
        '202S-1_ACCTGA_L001_R1_001.fastq.gz'
    and you want to split it into several output files with my format, use:
        's_G1_L001_R1_%03d.fastq.1.gz'
    '''
    kwargs['tempdir'] = output_dir
    print "Saving temp files to " + kwargs['tempdir']
    for chunk_num, chunk in enumerate(chunk_file(input_fn, **kwargs)):
        print "\tMoving file from %s to %s" % (chunk, output_form % chunk_num)
        shutil.move(chunk, output_form % chunk_num)
    
