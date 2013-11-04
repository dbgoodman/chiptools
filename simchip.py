import random
import itertools
import string
import warnings
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import fasta

BIN_PATH = '/Users/dbgoodman/Dropbox/coli regulatory elements/code/bin'
FASTA_FN = '/Users/dbgoodman/Dropbox/coli regulatory elements/'+\
           'library design/SimpleMatrix.fasta'
#          'library design/FinalLibraryAgilentUpload.fasta'
#          'library design/FinalCCDSShuffle.fasta'
           
LIBRARY_SIZE = 200000 #200 thousand colonies

###
# The Custom Read 1 adapter includes the NdeI cut site so we will not read 
# those bases, but we since we are using the default adapter for the AscI site
# we want to keep those bases on the fragment. 
#
# We want to keep the restriction sequence on the AscI side (read 2):
#   CATATG
# and remove the restriction sequence on the NdeI side (read 1):
#   GGCGCGCC
#
#SEQ UP TO NOT INCLUDING FIRST RESTRICTION SITE (AscI, GGCGCGCC)
STA_REMOVE = 'AATCCTTGCGTCAATGGTTC'  #for SimpleMatrix
    #OR TGTCGTGCCTCTTTATCTGT for CDSShuf
#SEQ AFTER AND INCLUDING SECOND RESTRICTION SITE (NdeI, CATATG)
END_REMOVE = 'CATATGCGTGTAAAATCCGAGAACCC' #for SimpleMatrix
    #OR revcom(GCTTCGGTGTATCGGAAATGCATATG) for CDSShuf

###
# Because the reads are in the opposite orientation we want to take the
# reverse complement after we trim the primer sequence off the reads.
#USE_REVCOM = True <- Automatic, no flag required

###
# append a constant sequence and a 0-5 bp variable sequence to the end of
# each fragment:
FRAG_END_CONST = 'ATGACTAAGCTTTTCATTGTC'
FRAG_END_VAR = ['', 'A', 'TA', 'ATG', 'ATGC', 'AATGC']
FRAG_END_APPEND = lambda: ''.join([FRAG_END_CONST,
                                   random.choice(FRAG_END_VAR)])
                                   
PCT_TRV = .0036
PCT_TRI = PCT_TRV+.0013
PCT_DEL = PCT_TRI+.0011
PCT_MATCH = 1
TRANSITION_TABLE = string.maketrans('ATGC','GCAT')
    
class mutate_seq:
    def __init__(self):
        #probabilities for transversion, transition, deletion, match
        self.prob_list = [[mutate_seq.do_transversion, PCT_TRV],
                          [mutate_seq.do_transition, PCT_TRI],
                          [mutate_seq.do_deletion, PCT_DEL],
                          [mutate_seq.do_nothing, PCT_MATCH]]
        
        #convert probs to be along a uniform distribution, 0 to 1
        start_dist = 0
        
        for i, (fxn, pct) in enumerate(self.prob_list):
            self.prob_list[i][1] = start_dist
            start_dist = start_dist + pct
    
    @staticmethod
    def do_transversion(letter):
        if letter == 'A' or letter == 'G':
            return random.choice(['C','T'])
        elif letter == 'C' or letter == 'T':
            return random.choice(['A','G'])
    
    @staticmethod
    def do_transition(letter):
       return letter.translate(TRANSITION_TABLE)
    
    @staticmethod
    def do_deletion(letter):
        return ''
    
    @staticmethod
    def do_nothing(letter):
        return letter
    
    def mutate_letter(self, letter):
        rand = random.random()
        i = 0
        while i < 3 and rand >= zip(*self.prob_list)[1][i]:
            i += 1
        return self.prob_list[i][0](letter)
    
    def mutate_sequence(self, seq):
        return Seq(''.join([self.mutate_letter(i) for i in seq]))

def trim_and_append(fasta_seq):
    '''
    simulate the restriction digest and append illumina library sequence
    '''
    trim_slice = slice(len(STA_REMOVE),-len(END_REMOVE))
    fasta_revcom = fasta_seq[trim_slice].reverse_complement()
    return ''.join([fasta_revcom.tostring(),FRAG_END_APPEND()])
    
def gen_clone_library(fasta_unif_dist, num_clones= LIBRARY_SIZE, **kwargs):
    '''
    generate a library of clones (i.e. cells) with mutated sequences. These
    will be mutated just like chip oligos.
    '''
    
    clone_library = []
    mutator = mutate_seq()
    
    clone_library = []
    genned_fastas = (fasta_unif_dist.next() for i in range(num_clones))
    for i, name_seq_pair in enumerate(itertools.izip(*[iter(genned_fastas)])):
        name,seq = name_seq_pair[0]
        trimmed_sr = fasta.extract_between_restr_sites(
                                seq.upper(),
                                STA_REMOVE,
                                END_REMOVE)
        if not trimmed_sr: 
            raise(InputError,"Library member does not specified trim sites!")
            
        trimmed_sr.seq = mutator.mutate_sequence(trimmed_sr.upper())
        trimmed_sr.seq = Seq(''.join([trimmed_sr.seq.tostring(),
                                  FRAG_END_APPEND()]))

        record = SeqRecord(trimmed_sr.seq.reverse_complement(),
                           id="fragment_"+str(i),
                           description=name)
                           
        clone_library.append(record)
    return clone_library

def sim_pcr_amp(clone_library, num_fragments= 5000000, output= sys.stdout, 
    **kwargs):
    '''
    simulates PCR / clonal growth of the library by randomly picking 
    mutated members from the library until all fragments have been generated
    '''
    lib_size = len(clone_library)
    for i in range(num_fragments):
        j = random.randint(0,lib_size-1)
        print >> output, clone_library[j].format('fasta')

def simulate_chip(chip_fasta_fn, num_fragments= 5000000, **kwargs):
    '''
    1. generate a mutated random library based on the fasta input
    2. amplify and send to simNGS
    '''
    #load fasta list
    fasta_list = fasta.load_fasta_file(chip_fasta_fn)
    #instatiate generator
    fasta_unif_dist = fasta.gen_seq_unif(fasta_list)
    print >> sys.stderr, 'Generating Clone Library...'
    #mutate sequences picked from generator as if they had synthesis errors
    clone_library = gen_clone_library(fasta_unif_dist, **kwargs)
    print >> sys.stderr, 'Generating Amplified Library...'
    #simulate clonal growth & PCR amplification of the mutated sequences
    sim_pcr_amp(clone_library, num_fragments, **kwargs)
    
if __name__ == "__main__":
    main()

