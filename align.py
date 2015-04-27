import subprocess
import itertools
import tempfile
import string
import os
import re
import sys

from collections import defaultdict
from StringIO import StringIO
from collections import defaultdict

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from cogent.app import cd_hit

import sample
import fasta

TRANSITION_TABLE = string.maketrans('ATGC','GCAT')
# ============================================================================
# = CD-HIT
# ============================================================================

def call_cdhit2d(library_file, query_file, output_file,
                 option_string='-p 1 -d 0 -g 1 -c 0.85 -S 10 -n 5 -t 5',
                 **kwargs):
    #cd-hit-2d -i library.fasta -i2 202-1.sanger_extracted.fasta
    #   -p 1 -d 0 -g 1 -c 0.85 -S 10 -n 5 -o 202-1.cd-hit-2d
    command = ' '.join([
        'cd-hit-2d ',
        '%s -i %s -i2 %s -o %s' % (option_string,
            library_file, query_file, output_file)])
    subprocess.call(command.split())
    print >> sys.stderr, 'Cluster file output to %s.clstr' % (output_file)

def cdhit2d(lib_seqrecords, query_seqrecords, output_prefix, **kwargs):

    #parse library directory
    if isinstance(lib_seqrecords,file):
        lib_seqrecords = fasta.extract_library(
            input_filename= lib_seqrecords,
            output=None, **kwargs)
    if not isinstance(lib_seqrecords[0],SeqRecord):
        raise(ValueError(
            'The library must be a FASTA file or a SeqRecord object!'))

    #parse query directory
    if isinstance(query_seqrecords,str):
        query_seqrecords = fasta.extract_sanger(
            query_seqrecords, output=None, **kwargs)

    if not isinstance(query_seqrecords[0],SeqRecord):
        raise(ValueError(
            'The query must be a directory or a SeqRecord object!'))

    lib_file = tempfile.NamedTemporaryFile(delete=False)
    SeqIO.write(lib_seqrecords, lib_file, 'fasta')
    lib_file.close()

    query_file = tempfile.NamedTemporaryFile(delete=False)
    SeqIO.write(query_seqrecords, query_file, 'fasta')
    query_file.close()

    #print 'LIB:',lib_file.name
    #print 'QRY:',query_file.name

    #call cdhit2d and output to 'output_prefix.clstr'
    call_cdhit2d(lib_file.name, query_file.name, output_prefix, **kwargs)

    #get cluster information
    (wells, lib_mems) = parse_cdhit_clusters(output_prefix+'.clstr',
        cull_empty= True, **kwargs)
    # lm_aligns = align_cdhit_pairs(lib_mems, lib_seqrecords,
    #                               query_seqrecords, **kwargs)

    cluster_stats(output_prefix, lib_file, query_file)

def align_cdhit_pairs(lib_clusters, librecords, seqrecords, **kwargs):
    print >> sys.stderr, 'Performing pairwise alignments...'
    #make lookup dictionaries from lists
    libdict = seqrecord_list_to_seqrecord_dict(librecords)
    seqdict = seqrecord_list_to_seqrecord_dict(seqrecords)
    aggregate_stats = defaultdict(int)
    total_pct_match = float(0)
    seq_count = 0
    lm_aligns = {}

    for lm,seqs in lib_clusters.items():
        lm_aligns[lm] = {}
        for seq in seqs:
            seq_count += 1

            if isinstance(seq,str):
                seq_label = seq
            else:
                seq_label = seq.id #for SangerPos namedtuples

            lm_aligns[lm][seq] = PairwiseAlign(libdict[lm],
                                                seqdict[seq_label], **kwargs)
    return lm_aligns

def parse_cdhit_clusters(cfname, output= None, cull_empty=False,
    plate_labels= True, **kwargs):

    print >> sys.stderr, 'Parsing Cluster file...'
    if isinstance(cfname,file):
        cfhandle = cfname
    else:
        cfhandle = open(os.path.join(cfname)).readlines()

    clines = [line.rstrip() for line in cfhandle]
    clstrs = cd_hit.parse_cdhit_clstr_file(clines)
    #cull clusters with one member
    if cull_empty: clstrs = filter(lambda i: len(i) > 1, clstrs)
    #print 2 formats: Seq\tWell Name\tLibrary Member
    #                 Library Member\tWell1\tWell2\tetc...

    #if well member follows regex (\d+-*\d*)-(\d+), it is a well
    queries = {}
    lib_mems = {}

    for clstr in clstrs:

        #NOTE: this assumes that second DB will always come first,
        #unsure if CD-HIT does this cannonically...(IT DOESNT)

        lib_mem = clstr.pop(0)
        lib_mems[lib_mem] = []
        for member in clstr:

            if plate_labels:
                sanger_pos = sample.parse_sanger_id(member)
                queries[sanger_pos] = lib_mem
                lib_mems[lib_mem].append(sanger_pos)
            else:
                queries[member] = lib_mem
                lib_mems[lib_mem].append(member)

    return queries, lib_mems

def read_stats(wells, lib_mems, lm_aligns=None, output= sys.stdout):

    #these will stay blank if lm_aligns is not present
    count_cols = ''
    dcount_cols = ''
    count_str = ''
    dcount_str = ''
    pct_match = ''
    header_printed = False

    #print header
    header = ['num','well','seq','% match']

    for well in sorted(wells,key=lambda well: int(well.num) ):

        lib_mem = wells[well]

        if lm_aligns:
            pwa = lm_aligns[lib_mem][well]

            #set up count_cols if this is the first line
            if count_cols == '':
                count_cols = PairwiseAlign.count_cols
                dcount_cols = PairwiseAlign.detailcount_cols

            #make the pct match string
            pct_match = '%3.1f%%' % (lm_aligns[lib_mem][well].pct_match * 100)

            #make the count string
            count_str = '\t'.join([str(pwa.count[col])
                                   for col in count_cols])
            dcount_str = '\t'.join([str(pwa.detailcount[col])
                                    for col in dcount_cols])

        #if we haven't printed the header yet:
        if not header_printed:
            print >> output, '\t'.join(header+count_cols+dcount_cols)
            header_printed = True

        #finally, print the row
        print >> output, '\t'.join([
            well.id, well.well, lib_mem, pct_match, count_str, dcount_str])

def lib_stats(lib_mems, output= sys.stdout):
    #print matches by library member:
    for lib_mem in sorted(lib_mems):
        print >> output, '\t'.join(
            [lib_mem]+[m.id for m in lib_mems[lib_mem]])

def align_stats(lm_aligns,librecords,seqrecords, output= sys.stdout):
    #make lookup dictionaries from lists
    libdict = seqrecord_list_to_seqrecord_dict(librecords)
    seqdict = seqrecord_list_to_seqrecord_dict(seqrecords)

    for lma, seqs in lm_aligns.items():
        for seq, align_obj in seqs.items():
            print >> output, align_obj

def cluster_stats(output_prefix, cluster_file, lib_seq, query_seq, **kwargs):

    #parse library
    if isinstance(lib_seq,str) or isinstance(lib_seq,file):
        lib_seqrecords = fasta.load_fasta_file(lib_seq)
    elif isinstance(lib_seq[0],SeqRecord):
        lib_seqrecords = lib_seq
    else:
        raise(ValueError(
            'The library must be a FASTA file or a SeqRecord object!'))

    #parse queries
    if isinstance(query_seq,str) or isinstance(query_seq,file):
        query_seqrecords = fasta.load_fasta_file(query_seq)
    elif isinstance(query_seq[0],SeqRecord):
        query_seqrecords = query_seq
    else:
        raise(ValueError(
            'The query must be a FASTA file or a SeqRecord object!'))

    #get cluster information
    #print kwargs
    (queries, lib_mems) = parse_cdhit_clusters(cluster_file, **kwargs)
    lm_aligns = align_cdhit_pairs(lib_mems, lib_seqrecords,
                                  query_seqrecords, **kwargs)

    #print statistics to files

    if kwargs['plate_labels']:
        cstats_file = open(output_prefix+'.cstats', 'w')
        read_stats(queries, lib_mems,
            lm_aligns= lm_aligns, output= cstats_file)

    libstats_file = open(output_prefix+'.libstats', 'w')
    lib_stats(lib_mems, output=libstats_file)

    aligns_file = open(output_prefix+'.aligns', 'w')
    align_stats(lm_aligns,lib_seqrecords,query_seqrecords, output=aligns_file)

def seqrecord_list_to_seqrecord_dict(seqrecord_list):
    return dict([(sr.id,sr) for sr in seqrecord_list])

# ============================================================================
# = PAIRWISE ALIGN
# ============================================================================

class PairwiseAlign:

    def __init__(self, seq1, seq2, local= False, from_multi= False,
                is_coding=0, **kwargs):

        #print local
        #print kwargs

        #deal with multiple object types:
        if isinstance(seq1,SeqRecord):
            self.id1 = seq1.id
            self.sr1 = seq1
            self.seq1 = self.sr1.seq
        elif isinstance(seq1,Seq):
            self.id1 = 'Sequence1'
            self.sr1 = None
            self.seq1 = seq1
        elif isinstance(seq1,str):
            self.id1 = 'Sequence1'
            self.sr1 = None
            self.seq1 = Seq(seq1)
        elif isinstance(seq1,tuple) and len(seq1) == 2:
            self.id1 = seq1[0]
            self.sr1 = None
            self.seq1 = Seq(seq1[1])

        if isinstance(seq2,SeqRecord):
            self.id2 = seq2.id
            self.sr2 = seq2
            self.seq2 = self.sr2.seq
        elif isinstance(seq2,Seq):
            self.id2 = 'Sequence2'
            self.sr2 = None
            self.seq2 = seq2
        elif isinstance(seq2,tuple) and len(seq2) == 2:
             self.id2 = seq2[0]
             self.sr2 = None
             self.seq2 = Seq(seq2[1])

        #if this is from a multiple alignment (via MUSCLE) then get the
        #aligned sequences here. Otherwise, do the all-python slow alignment.
        if from_multi:
            self.align1 = kwargs['align1']
            self.align2 = kwargs['align2']
            self.begin = None
            self.end = None
            self.score = None

        else:
            #get align type
            if local:
                aligntype = pairwise2.align.localdd
            else:
                aligntype = pairwise2.align.globaldd

            (self.align1, self.align2,
             score, begin, end) = aligntype(
                self.seq1.upper().tostring(),
                self.seq2.upper().tostring(),
                self.match_dict(),
                -5, -.5, -1, -.1, one_alignment_only=True)[0]

        self.count = defaultdict(int)
        self.detailcount = defaultdict(int)

        #add up matches and mismatches
        for i,j in zip(self.align1,self.align2):
            if i == j:
                if i != '-':
                    self.count['match'] += 1
                    self.detailcount['match_'+i] += 1
            elif j == 'N':
                self.count['match'] += 1
                self.detailcount['match_N'] += 1
            elif i == '-':
                self.count['ins'] += 1
                self.detailcount['ins_'+j] += 1
            elif j == '-':
                self.count['del'] += 1
                self.detailcount['del_'+i] += 1
            elif i.translate(TRANSITION_TABLE) == j:
                self.count['tsit'] += 1
            else:
                self.count['tvers'] += 1

        self.pct_match = (float(self.count['match'])
            / float(max(len(self.seq1),len(self.seq2))))

        self.is_coding = is_coding

        self.lev = self.levenshtein()

        if self.is_coding != 0:
            self.translate_align()


    @staticmethod
    def match_dict():
        pairs = list(itertools.product('ATGCN', repeat=2))
        match_dict = {}

        for p in pairs:
            #matches worth 1
            if p[0] == p[1]:
                match_dict[p] = 1
            #Ns count as match
            elif p[0] == 'N' or p[1] == 'N':
                match_dict[p] = 1
            #mismatches count as -2
            elif p[0] != p[1]:
                match_dict[p] = -2

        return match_dict

    @staticmethod
    def match_cmp_str(align1,align2):
        match_cmp_str = []
        for i,j in zip(align1, align2):
            if i != j: match_cmp_str.append('*')
            else: match_cmp_str.append(' ')
        return ''.join(match_cmp_str)

    def translate_align(self):
        start_pos = self.is_coding - 1

        self.trans1 =  self.seq1[start_pos:].translate(
                                                table="Bacterial").tostring()
        self.trans2 =  self.seq2[start_pos:].translate(
                                                table="Bacterial").tostring()

        untranslated =  ' ' * (start_pos)
        self.spaced_tr1 = untranslated+'  '.join([ aa for aa in self.trans1])
        self.spaced_tr2 = untranslated+'  '.join([ aa for aa in self.trans2])

    def __str__(self):
        output = StringIO()

        print >> output, "%-28s%s" % (self.id1, self.align1)
        print >> output, "%-28s%s" % (' ',
            PairwiseAlign.match_cmp_str(self.align1,self.align2))
        print >> output, "%-28s%s" % (self.id2, self.align2)

        if self.is_coding:
            print >> output, "%-28s%s" % (self.id1, self.spaced_tr1)
            print >> output, "%-28s%s" % (' ',
                PairwiseAlign.match_cmp_str(self.spaced_tr1,self.spaced_tr2))
            print >> output, "%-38s%s" % (self.id2, self.spaced_tr2)

        print >> output, ''.join(['%3.1f %%' % (self.pct_match*100)]
                                 + [' %s:%d' % (k,v)
                                    for k,v in self.count.items()]
                                 + [' %s:%d' % (k,v)
                                    for k,v in self.detailcount.items()])
        print >> output, '-'*(20+len(self.align1)+5)

        return output.getvalue()

    def levenshtein(self):
        ''' this function computes the edit distance between the two strings,
            but this should only need to happen once. it is stored in the lev
            member variable upon __init__()
        '''

        s1 = self.seq1.upper()
        s2 = self.seq2.upper()

        if len(s1) < len(s2):
            s2 = self.seq1.upper()
            s1 = self.seq2.upper()
        if not s1:
            return len(s2)

        previous_row = xrange(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
                deletions = current_row[j] + 1       # than s2
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row

        return previous_row[-1]

    #define count names here:
    count_cols = ['match','ins','del','tsit','tvers']
    detailcount_cols = (['match_'+n for n in 'ATGCN']
                        + ['del_'+n for n in 'ATGC'])



