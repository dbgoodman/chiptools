import gzip
import tempfile
import subprocess
import sys
import os
import warnings
import re
import pprint
import tempfile
import shutil
import glob
import time
import math
import itertools
import pickle
import shutil

from multiprocessing import Pool
from collections import defaultdict, namedtuple
from StringIO import StringIO
from shove import Shove

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC

import fasta
import counter
import align
import util
import pll_reduce


##STEPS########
# ==========run seqprep:=========
# - output files for r1, r2, failed and merged reads
# - fwd and reverse adapter sequences
# - X 1

# SeqPrep -f 202_R1.fq.gz -r 202_R2.fq.gz \
#   -1 $seqdir/202_R1.sp.fq.gz \
#   -2 $seqdir/202_R2.sp.fq.gz  \
#   -3 $seqdir/202_R1.sp.d.fq.gz \
#   -4 $seqdir/202_R2.sp.d.fq.gz  \
#   -s $seqdir/202_M.sp.fq.gz \
#   -A GGCGCGCCATGACTAAGCTTTTCATTGTCATGC \
#   -B CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT \
#   -X 1
# ==========create fasta with counts=========
# from the seqprep output, count fragments per sequence (with a hash) (in python)
# ==========cluster against library=========
# with another FASTA file (library), compare reads against the library and find
# the closest library member (using cd-2-hit) 
# ... get an edit distance for each using levenshtein distance (in python)

# a list of differences (insertions, mismatches, deletions) similar to CIGAR

# give the seq a name: original_seq.edit_distance.arbitrary_number

def call_seqprep(output_prefix, fq_read1_fn, fq_read2_fn,
    adapter1, adapter2, **kwargs):
    
    if not kwargs['multi']:
        print >> sys.stderr, 'Merging FASTQ reads with SeqPrep...'
        sys.stderr.flush()
    
    #generate file names
    files = {
    'fq1_out_fn': output_prefix+'.R1.sp.fq.gz',
    'fq2_out_fn': output_prefix+'.R2.sp.fq.gz',
    'fq1_disc_fn': output_prefix+'.R1.disc.fq.gz',
    'fq2_disc_fn': output_prefix+'.R2.disc.fq.gz',
    'merged_out_fn': output_prefix+'.M.fq.gz'}
    
    #skip if the option is turned on.
    if kwargs['skip_finished'] and os.path.exists(files['merged_out_fn']):
        if not kwargs['multi']:
            print >> sys.stderr, '\tFiles present; SeqPrep skipped.\n'
        return files
    
    #make command 
    command = '''
        SeqPrep 
        -f %(fq_read1_fn)s
        -r %(fq_read2_fn)s
        -1 %(fq1_out_fn)s
        -2 %(fq2_out_fn)s
        -3 %(fq1_disc_fn)s
        -4 %(fq2_disc_fn)s
        -s %(merged_out_fn)s
        -A %(adapter1)s -B %(adapter2)s
        -X 1''' % (dict(locals().items() + files.items()))
    
    #split command and wait for it to finish, pipe stderr and stdout to obj
    child = subprocess.Popen(command.split(),
        stdout= subprocess.PIPE,
        stderr= subprocess.PIPE)
    child.wait()
    
    #get stats, add to file dict for parsing later
    stats = {}
    for line in child.stderr.readlines():
        matches = re.match(r'(Pairs[ \w]+|CPU)[\w() ]*:\t([\d.]+)',line)
        if matches:
            stats[matches.groups()[0]] = util.num(matches.groups()[1])
    files['stats'] = stats
    
    if not kwargs['multi']:
        print >> sys.stderr, '\tSeqPrep output to %s.*\n' % (output_prefix)
    return files

def ngs_seqprep(output_prefix, fq_read1_fn, fq_read2_fn, adapter1, adapter2,
    **kwargs):
    
    files = call_seqprep(output_prefix, fq_read1_fn, fq_read2_fn,
        adapter1, adapter2, **kwargs)
    
    if kwargs['no_manual_trim']:
        return files
    
    files['trimmed_out_fn'] = output_prefix+'.MT.seq.gz'
    
    #skip if the option is turned on.
    if kwargs['skip_finished'] and os.path.exists(files['trimmed_out_fn']):
        if not kwargs['multi']:
            print >> sys.stderr, '\tFiles present; manual trim skipped.\n'
        return files
    
    #cull remaining sequences whose adapters were not trimmed
    trim_cmd = ' '.join([
        'cat <(zcat %(merged_out_fn)s | grep -P \'^[ATGC]+$\'',
             '| grep -Pv \'(?<=%(rs1)s)[ATGC]+(?=%(rs2)s)\')',
             '<(zcat %(merged_out_fn)s | grep -P \'^[ATGC]+$\'',
             '| grep -Po \'(?<=%(rs1)s)[ATGC]+(?=%(rs2)s)\')',
             '| sort | gzip > %(trimmed_out_fn)s']) % (
                dict(kwargs.items() + files.items()))

    subprocess.call(trim_cmd, shell=True, executable='/bin/bash')
    
    assert len(files) > 0
    
    return files

# ============================================================================
# = SINGLE FASTQ INPUT MODE
# ============================================================================

def load_single_input_files(fq_read1, fq_read2, lib_records, **kwargs):
    '''load/extract fastq/fasta files:
        we want to potentially trim the library at cut sites, but leave
        the reads alone (since bringing them into python would be slow)
        we just want to check that the fastq files look right.'''
    
    #parse library file/object
    if isinstance(lib_records,file):
        lib_records_fn = lib_records.name
        lib_records = fasta.extract_library(
            input_filename= lib_records, 
            output=None, **kwargs)
        
    if not isinstance(lib_records[0],SeqRecord):
        raise(ValueError(
            'The library must be a FASTA file or a SeqRecord object!'))                
        
    #check fastq files
    for fastq in fq_read1, fq_read2:
        if not isinstance(fastq,file) and fastq.endswith('.gz'):
            warnings.warn(fastq+
            ' must be a gzipped FASTQ file but doesn\'t end in .gz...')
    
    print >> sys.stderr, '\tSequence Extraction complete.\n'
    sys.stderr.flush()
    
    return lib_records

# ============================================================================
# = MULTIPLE FASTQ INPUT MODE
# ============================================================================

def ngs_multi_input_files(fq_reads, regex_str, lib_records, 
    adapter1, adapter2, output_prefix, **kwargs):
    
    '''load/extract fastq/fasta files:
        we want to potentially trim the library at cut sites, but leave
        the reads alone (since bringing them into python would be slow)
        we just want to check that the fastq files look right.
        
        We also want to group the reads by tile, read number, and bin.'''
        
    
    #fq_tree[bin][tile][read]
    fq_tree = defaultdict(lambda: defaultdict(dict))
    
    #create a tree of all the fastq files
    for fq_file in fq_reads:
        fq_components = re.match(regex_str, fq_file.name)
        fq_read = fq_components.group('read')
        fq_bin = fq_components.group('bin')
        fq_tile = fq_components.group('tile')
        
        fq_tree[fq_bin][fq_tile][fq_read] = fq_file
    
    #check the tree. print out the number of bins, the number of tiles, and
    #make sure that each bin has the same tiles and each tile
    #has only two reads, 1 and 2. 
    
    #check that tiles are the same
    all_tiles = [frozenset(bin.keys()) for bin in fq_tree.values()]
    if not len(set(all_tiles)) <= 1:
        raise( ValueError('Tiles are not the same among bins: '
            + print_fq_tree(fq_tree)))
    
    #check that each tile has only two reads, called 1 and 2:     
    all_reads = [
        frozenset([frozenset(bin[tile].keys()) for tile in bin.keys()])
            for bin in fq_tree.values()]
    if not len(set(all_reads)) <= 1:
        raise( ValueError('Not all tiles have the same two read files: '
            + print_fq_tree(fq_tree)))
            
    if not set(all_reads).pop() == frozenset([frozenset(['1', '2'])]):
        raise( ValueError('All reads must be named 1 and 2: '
            + print_fq_tree(fq_tree)))
    
    #parse library file/object
    if isinstance(lib_records,file):
        lib_records_fn = lib_records.name
        lib_records = fasta.extract_library(
            input_filename= lib_records, 
            output=None, **kwargs)
        
    if not isinstance(lib_records[0],SeqRecord):
        raise(ValueError(
            'The library must be a FASTA file or a SeqRecord object!'))
        
    print >> sys.stderr, 'Sending SeqPrep jobs to workers...'
    sys.stderr.flush()

    pool = Pool(None)
    sp_jobs = []
    
    for bin in fq_tree.keys():
        
        for tile, files in fq_tree[bin].items():
            
            files['sp_files'] = {}
            
            fq_read1_fn = files['1'].name
            fq_read2_fn = files['2'].name
                                    
            #tile and bin prefix
            t_b_prefix = output_prefix+'.%s.%s' % (bin, tile)
            
            sp_jobs.append(pool.apply_async(ngs_seqprep,
                (t_b_prefix, fq_read1_fn, fq_read2_fn, 
                adapter1, adapter2), kwargs, 
                callback= files['sp_files'].update))
            
    #after submitting all tiles, wait for them all to finish. 
    job_count = 0
    while True:
        time.sleep(1)
        completed = sum([spj.ready() for spj in sp_jobs])
        if completed > job_count:
            print >> sys.stderr, '\t%d / %d SeqPrep jobs completed.' % (
                completed, len(sp_jobs))
            sys.stderr.flush()
            job_count = completed

        if all([spj.ready() for spj in sp_jobs]):
            break
            
    #wait and close old pool for seqprep, start new pool
    pool.close()
    pool.join()
    print >> sys.stderr, '\tSeqprep complete.\n'
    sys.stderr.flush()
    
    
    #print out the sp_stats file to get number of reads trimmed, discarded,
    # etc. If we are skipping finished parts of the pipeline, then we will not
    #recompute the seqprep stats either...
    if not kwargs['skip_finished']: sp_stats(fq_tree, output_prefix)
    
    return fq_tree, lib_records

def ngs_count_unique_reads(seq_file, **kwargs):
    ''' for counting, use the shell to count unique reads per
        tile/bin, then return a file-like obj with the counts for that bin.'''
    
    seq_file = util.open_zcat(seq_file)
    
    count_file = tempfile.NamedTemporaryFile(delete=False)
    
    command = 'grep -Pi \'^[ANTGC]+$\''
    
    if kwargs['revcom']:
        command += ' | perl -ne \'chomp; $_ =~ tr/ACGTacgt/TGCAtgca/;'
        command += ' print reverse($_)."\n";\''
    
    command += ' | sort | uniq -c | gzip'
    
    subprocess.call(command,
        stdout=count_file,
        stdin=seq_file,
        shell=True,
        executable='/bin/bash') 
    
    return count_file.name

def ngs_multi_count_reads(fq_tree, output_prefix, **kwargs):
    
    print >> sys.stderr, 'Sending count jobs to workers...'
    sys.stderr.flush()
    
    counts_filename = output_prefix+'.shove.dict'
    
    if kwargs['skip_finished'] and os.path.exists(counts_filename):
        print >> sys.stderr, '\tFile present; counting skipped.\n'
        sys.stderr.flush()
        return
    
    pool = Pool() #open a worker pool
    bin_tile_files = {} #seqprep read files per bin 
    bin_tile_files['all'] = [] #(and total under 'all')
    bin_count_files = {} #resulting files per bin 
    count_jobs = {} #keep track of parallel jobs
    
    #use a different file name if we've manually trimmed the sequences    
    no_trim = kwargs['no_manual_trim']
    file_key = 'merged_out_fn' if no_trim else 'trimmed_out_fn'
    
    #get the read files
    for bin in fq_tree.keys():
        
        tile_files = [fq_tree[bin][tile]['sp_files'][file_key] 
            for tile in fq_tree[bin].keys()]
        
        bin_tile_files[bin] = tile_files
        bin_tile_files['all'].extend(tile_files)
    
    #count reads in total
    
    #print bin_tile_files['all']
    
    count_jobs['all'] = pool.apply_async(ngs_count_unique_reads,
        [bin_tile_files['all']], kwargs)
    
    #count reads per bin
    for bin in fq_tree.keys():
        tile_files = [fq_tree[bin][tile]['sp_files'][file_key] 
            for tile in fq_tree[bin].keys()]
            
        count_jobs[bin] = pool.apply_async(ngs_count_unique_reads, 
            [tile_files], kwargs)
    
    #after submitting all tiles, wait for them all to finish. 
    job_count = 1
    while True:
        time.sleep(1)
        completed = sum([cj.ready() for cj in count_jobs.values()])
        if completed > (job_count):
            print >> sys.stderr, '\t%d / %d workers completed.' % (
                completed, len(count_jobs))
            sys.stderr.flush()
            job_count = completed

        if all([cj.ready() for cj in count_jobs.values()]): 
             for bin, cj in count_jobs.items():
                 bin_count_files[bin] = cj.get() 
             break

    pool.close()
    pool.join()
    print >> sys.stderr, '\tBin counting complete.\n'
    sys.stderr.flush()
    
    print >> sys.stderr, 'Gathering counts from bins (not parallelized)...'
    sys.stderr.flush()
    
    #now, we want to open file objects for all bins and create a dict of lists
    #where the dict keys are seqs and each value is a list for each bin and
    #'all' for all bins. We also need to assert that 'all' is the sum of every
    #bin.
    
    #open all filehandles
    bin_count_fhandles = {}
    #for holding current lines per file
    bin_count_clines = {}
    
    for bin, fn in bin_count_files.items():
        bin_count_fhandles[bin] = util.open_zcat(fn)
    
    #dict for sequences
    fq_counts = {}

    #default dict for bins and all
    empty_dict = dict(zip(bin_count_fhandles, [0] * len(bin_count_fhandles)))
        
    total_count_fhandle = bin_count_fhandles['all']
    del bin_count_fhandles['all'] #remove all so we can loop through the bins
    
    #SeqCount namedtuple
    SeqCount = namedtuple('SeqCount', ['bin', 'count', 'seq'])
    make_sc = lambda bin, line: SeqCount(bin, *line.split())
    
    
    #for every bin, grab the first line and put it in a 'current lines' dict
    for bin, fh in bin_count_fhandles.items():
        bin_count_clines[bin] = make_sc(bin, fh.next())
    
    for line in total_count_fhandle:
        count_total = 0 #for checking that 'all' == sum(bins)
        total_seqcount = make_sc('all', line) 
        fq_counts[total_seqcount.seq] = empty_dict.copy()
        fq_counts[total_seqcount.seq]['all'] = int(total_seqcount.count)
        
        for bin, seqcount in bin_count_clines.items():
            if total_seqcount.seq == seqcount.seq:
                fq_counts[total_seqcount.seq][bin] = int(seqcount.count)
                count_total += int(seqcount.count)  
            
                try:
                    bin_count_clines[bin] = make_sc(
                        bin, bin_count_fhandles[bin].next())
                except StopIteration:
                    del bin_count_clines[bin]
        
        
        assert fq_counts[total_seqcount.seq]['all'] == count_total
    
    #finally, pickle this fq_counts dict
    counts_file = open(counts_filename, 'wb')
    pickle.dump(fq_counts, counts_file)
    #fq_shove = Shove('file://'+counts_filename)
    #fq_shove.update(fq_counts)
    
    print >> sys.stderr, '\tCount object made successfully.\n'
    sys.stderr.flush()

def ngs_unique_records(output_prefix, bins=False, **kwargs):
    '''make a seqrecord array for each unique seq, store count(s) for each,
       write to a fasta file'''
       
    fqcounts_filename = output_prefix+'.shove.dict'
    fq_counts = pickle.load(open(fqcounts_filename, 'rb'))
    #fq_counts = Shove('file://'+fqcounts_filename)
    
    fq_unique_records = []    
    print >> sys.stderr, 'Writing records for each unique sequence...'
    sys.stderr.flush()
    
    counts_filename = output_prefix+'.counts.fa'
    
    if kwargs['skip_finished'] and os.path.exists(counts_filename):
        print >> sys.stderr, '\tFile present; record writing skipped.\n'
        sys.stderr.flush()
        return counts_filename

    counts_file = open(counts_filename, 'w')
    
    for member_id, member_tuple in enumerate(fq_counts.items()):
        unique_member, member_count = member_tuple

        member_annotation = {
            'member_id': member_id,
            'matched': False,
            'designed': None}
    
        #if bins is true, then each value in the fq_counts will be a dict of 
        #bins with an 'all' key corresponding to total counts
        if bins:
            bin_counts = member_count
            member_annotation['count'] = bin_counts['all']
            del bin_counts['all']
        
            for bin, count in bin_counts.items():
                member_annotation['bin.'+bin] = count
        else:
            member_annotation['count'] = member_count
    
        if kwargs['multi']:
            description='unmatched'+' '+member_annotation.__repr__()
            seq_id='U'+str(member_id)
            print >> counts_file, ('>%(seq_id)s %(description)s\n'
                +'%(unique_member)s') % dict(locals().items())
        else:
            fq_unique_records.append(
                SeqRecord(seq=Seq(unique_member,IUPAC.unambiguous_dna),
                    id='U'+str(member_id),
                description='unmatched'+' '+member_annotation.__repr__(),
                annotations=member_annotation))
        
    if not kwargs['multi']:
        SeqIO.write(fq_unique_records, counts_file, 'fasta')
    
    print >> sys.stderr, '\tRecord printing complete.\n'
    sys.stderr.flush()
    
    if kwargs['multi']:
        counts_file.close()
        return counts_filename
    else:
        return fq_unique_records


def ngs_cdhit2d(lib_records, output_prefix, **kwargs):
    '''run cdhit2d on whole unique read library against designed library'''
    
    print >> sys.stderr, (
        'Aligning unique sequences to designed library with cd-hit-2d...')
    sys.stderr.flush()
    
    counts_filename = output_prefix+'.counts.fa'
    
    if kwargs['skip_finished'] and os.path.exists(output_prefix+'.clstr'):
        print >> sys.stderr, '\tFiles present; cd-hit-2d skipped.\n'
        sys.stderr.flush()
    else:
        if kwargs['multi']:
            query_file = open(counts_filename, 'r')
            query_file.close()
        else:
            query_file = tempfile.NamedTemporaryFile(delete=False)
            SeqIO.write(fq_unique_records, query_file, 'fasta')
            query_file.close()

        lib_file = tempfile.NamedTemporaryFile(delete=False)
        SeqIO.write(lib_records, lib_file, 'fasta')
        lib_file.close()
        
        #option_string='-G 1 -aL .9 -p 1 -d 0 -g 1 -c 0.8 -aS 1 -n 5 -t 5 -T 0'
        option_string='-G 1 -aS 1 -aL 1 -d 0 -g 1 -c 1 -n 5 -t 5 -T 0 -M 7500'
        
        align.call_cdhit2d(lib_file.name, query_file.name, output_prefix, 
            option_string=option_string, **kwargs)
        
    (queries, lib_mems) = align.parse_cdhit_clusters(
        output_prefix+'.clstr', plate_labels= False, **kwargs)

    print >> sys.stderr, (
        'Reading sequence dictionary objects...')
    sys.stderr.flush()
            
    fq_fhandle = open(counts_filename, 'r')
    
    #lm_shove = Shove('file://'+output_prefix + '.shove.libmems')
    lm_shove = Shove()
    #seqdict = Shove('file://'+output_prefix + '.shove.seqdict')
    seqdict = Shove()
    #seqannotdict = Shove('file://'+output_prefix + '.shove.sadict')
    
    seqdict = util.quick_fasta(counts_filename, 
        get_annotations= True, seqdict= seqdict)
    
    libdict = align.seqrecord_list_to_seqrecord_dict(lib_records)
    
    lm_shove.update(lib_mems)
    
    #for lm in lib_mems:
    #    lm_shove[lm] = lib_mems[lm]
    
    for member in queries:
        seqdict[member][1]['matched'] = True
        seqdict[member][1]['designed'] = queries[member]

    print >> sys.stderr, (
        '\tSeqRecord objects complete.\n')
    sys.stderr.flush()
    
    return libdict, seqdict

def ngs_muscle_align(libdict, seqdict, output_prefix, **kwargs):
    #if this designed library member has any unique sequences that match to it
    #then perform a pairwise alignment for each and update the relevant info
    #for that unique sequence. 
    #Also update the count for the unique member. 
    
    lib_mems = Shove('file://'+output_prefix + '.shove.libmems')
    #seqdict = Shove('file://'+output_prefix + '.shove.seqdict')
    
    #====do alignment of clusters====
    
    #create parallel process pool with all threads
    pool = Pool()
    muscle_jobs = {}
    completed_jobs = 0
    align_kwargs ={'ref_idx': 0, 'local': False, 'is_coding': 0}
    print >> sys.stderr, 'Using MUSCLE to align clustered sequences...'
    sys.stderr.flush()    
    
    pairwise_aligns = Shove('file://'+output_prefix+'.shove.pwa')
    
    faalign_file = open(output_prefix+'.aligns.fa', 'w')
    
    def use_muscle_output(muscle_tuple):
        pwa, align_fa = muscle_tuple
        print >> faalign_file, align_fa
        pairwise_aligns.update(pwa)
    
    for lm, uniques in lib_mems.items():
        libdict[lm].annotations['count'] = len(uniques)
        
        #perform parallel multiple alignment with clustered seqs using MUSCLE                
        if len(uniques) == 0: continue
        
        muscle_jobs[lm] = pool.apply_async(align_recordlist, (
           [(lm,libdict[lm].seq.tostring(),)] 
           + [(uniq, seqdict[uniq][0],) for uniq in uniques],),
           {'ref_idx': 0, 'local': False, 'is_coding': 0})     
        
        # out_tuple = (align_recordlist(
        #     [(lm,libdict[lm].seqdict(),)] 
        #     + [(uniq, seqdict[uniq],) for uniq in uniques], **align_kwargs))
        # pairwise_aligns.update(out_tuple[0])
        # print >> faalign_file, out_tuple[1]
        #  
        # completed_jobs += 1
        # if completed_jobs % 100 == 0:
        #     print >> sys.stderr, '\t%d / %d clusters completed.' % (
        #         completed_jobs, len(lib_mems.items()))
        #     sys.stderr.flush()
    
    submitted_jobs = len(muscle_jobs.keys())
    print >> sys.stderr, '\t%d jobs submitted.' % (submitted_jobs)
    
    #wait for all jobs to finish        
    while True:
        for it in filter(lambda it: it[1].ready(), muscle_jobs.items()):
            lm, mj = it
            completed_jobs += 1
            if not mj.successful():
                import pdb
                pdb.set_trace()
            else:
                use_muscle_output(mj.get())
                del muscle_jobs[lm]
            if completed_jobs % 250 == 0:
                print >> sys.stderr, '\t%d / %d clusters completed.' % (
                    completed_jobs, submitted_jobs)
                sys.stderr.flush()
        
        if len(muscle_jobs) == 0:
            break
                        
    pool.close()
    pool.join()
        
    print >> sys.stderr, '\tMUSCLE alignments complete.\b'
    sys.stderr.flush()


# ============================================================================
# = PRINT STATS
# ============================================================================


#method for printing fq tree:
def print_fq_tree(fq_tree):
    ''' heirarchical tree of all fq.gz files inputted, used for debugging 
    purposes mostly.'''
    
    output = StringIO()
    
    for b in sorted(fq_tree.keys()):
        print >> output, 'Bin '+ b
        for t in sorted(fq_tree[b].keys()):
            print >> output, '  Tile ' + t
            for k,v in fq_tree[b][t].items():
                v_val = v.name if isinstance(v,file) else v
                print >> output, '     ' + k + ': ' + v_val
    
    return output.getvalue()

# ============================================================================

def sp_stats(fq_tree, output_prefix):
    ''' statistics file for each seqprep tile file, pairs merged, trimmed, 
    discarded, CPU times, etc.'''
    
    spstats_file = open(output_prefix+'.spstats', 'w')
    
    #print header
    print >> spstats_file, '\t'.join([
        'Bin',
        'Tile',
        'Processed',
        'Merged'
        'PctMerged',
        'Trimmed',
        'PctTrimmed',
        'Discarded',
        'PctDiscarded',
        'Time'])
    
    for bin in sorted(fq_tree.keys()):
        for tile,tdict in sorted(fq_tree[bin].items()):
            ppro = tdict['sp_files']['stats']['Pairs Processed']
            pmrg = tdict['sp_files']['stats']['Pairs Merged']
            pada = tdict['sp_files']['stats']['Pairs With Adapters']
            pdsc = tdict['sp_files']['stats']['Pairs Discarded']
            cput = tdict['sp_files']['stats']['CPU']
            
            #print line:
            print >> spstats_file, '\t'.join(map(str,[
                bin,tile,ppro,
                pmrg, pmrg/float(ppro),
                pada, pada/float(ppro),
                pdsc, pdsc/float(ppro),
                cput]))
    
    spstats_file.close()

# ============================================================================
    
def ngs_stats(seqdict, output_prefix, **kwargs):
    ''' statistics file for each unique read, for counts per bin, total counts
    differences from designed sequences, etc.'''
    
    #pairwise_aligns = Shove('file://'+output_prefix + '.shove.pwa')
    #seqdict = Shove('file://'+output_prefix + '.shove.seqdict')
    #seqannotdict = Shove('file://'+output_prefix + '.shove.sadict')
    
    print >> sys.stderr, 'Printing ustats file...'
    sys.stderr.flush()    
    
    ustats_file = open(output_prefix+'.ustats', 'w')
    
    # header = [
    #     "UniqRead",
    #     "Obs",
    #     "DsgnMatch",
    #     "LevDist",
    #     "#Dels",
    #     "#Ins",
    #     "#Mismatch",
    #     "PctMatch"]
    header = [
        "UniqRead",
        "Obs",
        "DsgnMatch"]
    
    #check and see if there are bins, and if so, add them to header
    bin_list = []
    if kwargs['multi']:
        for annot in seqdict.values()[0][1].keys():
            if annot.startswith('bin'):
                bin_list.append(annot)
        sorted_uniques = sorted(seqdict, 
            key=lambda u: seqdict[u][1]['count'], reverse= True)
    else:
        for annot in seqdict.values[0].annotations.keys():
            if annot.startswith('bin'):
                bin_list.append(annot)
        sorted_uniques = sorted(seqdict,
            key=lambda u: u.annotations['count'],
            reverse= True)
    
    bin_list.sort()
    header.extend(bin_list)
    
    #print header
    print >> ustats_file, '\t'.join(header)
    
    #print a line in the .ustats file for each unique sequence
    for uniq in sorted_uniques:
        
        #shortcuts to relevant variables
        uanno = seqdict[uniq][1]
        dlm = uanno['designed']
        if dlm == None:
            continue
        
        #designed library member
        #if uniq in pairwise_aligns:
            #ualign = pairwise_aligns[uniq]
            #uanno['matched'] = ualign.id1
            #uanno['align'] = ualign
            #uanno['lev'] = ualign.lev
            #uanno['pct_match'] = ualign.pct_match
            #uanno['al_counts'] = ualign.count
            #uanno['al_dcounts'] = ualign.detailcount            
            #dlm = uanno['designed']
            #lev = str(uanno['lev'])
            #pct = str(uanno['pct_match'])
            #dels = str(ualign.count['del'])
            #ins = str(ualign.count['ins'])
            #mm = str(ualign.count['tsit'] + ualign.count['tvers'])
        
        #else:
            #dlm = lev = pct = mm = dels = ins = 'NA'
        
        if bin_list:
            bin_vals = map(uanno.get,bin_list)
            bin_vals = map(str,bin_vals)
        else:
            bin_vals = []
        
        print >> ustats_file, '\t'.join([
            uniq,
            str(uanno['count']),
            dlm] + bin_vals)
        
        # print >> ustats_file, '\t'.join([
        #     uniq,
        #     str(uanno['count']),
        #     dlm,
        #     lev,
        #     dels,
        #     ins,
        #     mm,
        #     pct]
        #     + bin_vals)            
        
    ustats_file.close()

# ============================================================================
# = MAIN METHOD
# ============================================================================

def ngs_cluster(fq_read1, fq_read2, fq_files, regex_str,
    output_prefix, lib_records, adapter1, adapter2, **kwargs):
    
    #Single FASTQ Mode
    if fq_read1 and fq_read2 and not fq_files:
        print >> sys.stderr, 'NOTE: Working in Single FASTQ File Mode.\n'
        sys.stderr.flush()        
        
        kwargs['multi'] = False
        
        lib_records = load_single_input_files(fq_read1, fq_read2, lib_records)
        
        sp_files = ngs_seqprep(
            output_prefix, fq_read1.name, fq_read2.name, adapter1, adapter2, 
            **kwargs) 
        
        fq_counts = ngs_count_reads(
            sp_files['merged_out_fn'], output_prefix, **kwargs)
        
        fq_unique_records = ngs_unique_records(
            fq_counts, output_prefix, **kwargs)
        
    #Multi FASTQ Mode
    elif fq_files and not (fq_read1 and fq_read2):
        print >> sys.stderr, 'NOTE: Working in Multi FASTQ File Mode.\n'
        sys.stderr.flush()    
        if kwargs['skip_finished']:
            print >> sys.stderr, (
                'NOTE: Skipping steps where files are already present.\n')
            sys.stderr.flush()
        #else:
        #    for path in glob.glob( output_prefix + '*shove*' ): 
        #        try:
        #            shutil.rmtree( path )
        #        except:
        #            os.remove( path )
        
        kwargs['multi'] = True
        
        (fq_tree, lib_records) = ngs_multi_input_files(
            fq_files, regex_str, lib_records, adapter1, adapter2, 
            output_prefix, **kwargs)
        
        ngs_multi_count_reads(fq_tree, output_prefix, **kwargs)
                
        ngs_unique_records(output_prefix, bins=True, **kwargs)
        
    else:
        print >> sys.stderr, ('Please specify --read1 and --read2 files, or a'
                              +' list of FASTQ files with --fastq, not both.')
        
        sys.stderr.flush()
        sys.exit('Incorrect arguments.')
        
    libdict, seqdict = ngs_cdhit2d(lib_records, output_prefix, **kwargs)  
        
    #ngs_muscle_align(libdict, seqdict, output_prefix, **kwargs)
    
    #clear old temp files
    for filename in glob.glob( '/tmp/tmp*' ): os.remove( filename )
    
    ngs_stats(seqdict, output_prefix, **kwargs)


# ============================================================================
# = MULTIPLE ALIGN
# ============================================================================

def align_recordlist(clusterlist, ref_idx= 0, local= False, is_coding= 0, 
    **kwargs):
    '''
    Given a list of seqrecord objects, and specifying one as a reference 
    sequence, generate a list of pairwise align objects (see above) where the
    alignment is calculated from the multiple alignment over the entire set.
    This is much faster than using the all-python aligner above. 
    '''
    cline = MuscleCommandline(maxiters= 1, diags= True)
    child = subprocess.Popen(str(cline),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=(sys.platform!="win32"))

    for cluster_tuple in clusterlist:
        print >> child.stdin, '>%s\n%s' % (cluster_tuple)
    
    child.stdin.close()
    
    out1, out2 = itertools.tee(child.stdout)
    
    aligndict = util.quick_fasta(out1)
    
    (refid, refseq) = clusterlist[ref_idx]

    pairwise_aligns = {}

    #for each non-ref record, make a pairwise align record
    for idx,record in enumerate(clusterlist):
        if idx == ref_idx: continue
        seqid, seq = record
        
        #pairwise_aligns[seqid] = align.PairwiseAlign(
        #    clusterlist[ref_idx],
        #    record,
        #    align1 = aligndict[refid], 
        #    align2 = aligndict[seqid],
        #    from_multi= True, is_coding= is_coding, **kwargs)
    
    faa_output = ''.join(out2)
    
    child.stdout.close()
    
    return pairwise_aligns, faa_output

# ============================================================================
# = MAP REDUCE CODE, needs to be fixed:
# ============================================================================

# def ngs_multi_count_reads(fq_reads, regex_str, lib_records, 
#     adapter1, adapter2, output_prefix, **kwargs):
#     
#     #make a recursive pool so workers can have their own pools....
#     #http://stackoverflow.com/
#     # questions/6974695/python-process-pool-non-daemonic
#     outer_pool = pll_reduce.RecursivePool(8)
#     
#     mr_jobs = {} #list of map reduce jobs to hold async_results objects
#     
#     #make bin and total counters for reads
#     
#     no_trim = kwargs['no_manual_trim']
#     file_key = 'merged_out_fn' if no_trim else 'trimmed_out_fn'
#     for bin in fq_tree.keys(): 
#         tile_files = [fq_tree[bin][tile]['sp_files'][file_key] 
#             for tile in fq_tree[bin].keys()]
#                      
#         mr_jobs[bin] = outer_pool.apply_async(parallelize_bin_counting,
#             [bin, tile_files, output_prefix], kwargs)
#         #parallelize_bin_counting(bin, tile_files, output_prefix, **kwargs)
#     
#     bin_count_files = {}
#     while True:
#         if all([mrj.ready() for mrj in mr_jobs.values()]): 
#             [bin_count_files.update(mrj.get()) for mrj in mr_jobs.values()]
#             break
# 
#             
#     #we don't need this outer pool anymore...
#     outer_pool.close()
#     outer_pool.join()
#         
#     #once all bins are complete in this fashion, get total counts with another
#     #parallel reduce. First open all files for bin counts that we have copied
#     
#     #use another round of parallel reduce to get the count file for all bins
#     total_count_file = pll_reduce.parallel_reduce(
#         ngs_mcount_reduce, bin_count_files.values(), processes= 16)
#     
#     assert len(total_count_file) == 1
#     
#     total_count_file = total_count_file[0]
#     
#     #copy the file to its permanent location and delete the temp file 
#     total_count_filename = output_prefix +'.total.counts'
#     shutil.copy(total_count_file, total_count_filename)
#      
#     #We don't want to use ngs_unique_records any more because we don't want 
#     #counter objects for each bin; it would eat up too much memory.
#     #fq_unique_records = ngs_unique_records(total_counts, output_prefix, 
#     #    bins=bin_counts, **kwargs)
#     
#     #LEFT TO DO!
#     import pdb;
#     pdb.set_trace()
#     
#     return fq_unique_records
#     
#     
#     return lib_records, fq_unique_records, bin_counts

# def parallelize_bin_counting(bin, tile_filenames, output_prefix, **kwargs):
#     ''' This method allows us to us nested pools in our bin counting.
#         The outer pool is for different bins - this function is called
#         from that pool. I am using a parallelized map-reduce framework
#         to speed up sequence counting. First I bin per tile (the map
#         phase) then I combine the tiles by folding them together pairwise 
#         (the reduce phase). Both of these phases use the inner pools of
#         workers.'''
#         
#     bin_count_files = {} #will be merged with main dict on callback
#     
#     temp_count_filenames = []
# 
#     #map phase: generate sequence counts for individual tiles/bins
#     #==========
#     #non-recursive inner pool for map phase
#     tile_pool = Pool(4)
#     # a list to store worker async_result objects
#     tile_jobs = [] 
#         
#     for tfile in tile_filenames:
#     
#         #reads are named and stored differently if we manually trimmed 
#         #them after seqprep earlier            
#         tile_jobs.append(
#             tile_pool.apply_async(ngs_mcount_map,
#                 (tfile,),kwargs))
#         #temp_count_files.append(ngs_mcount_map(seq_file, **kwargs))
#         
#         #wait until all tile jobs complete
#         while True: 
#             if all([tj.ready() for tj in tile_jobs]): 
#                 temp_count_filenames = [tj.get() for tj in tile_jobs]
#                 break
#     
#     #close the inner pool for the tile map jobs
#     tile_pool.close()
#     tile_pool.join()
#         
#     #reduce: merge tiles per bin into bin counts
#     #=========
#     #inner pool will be made inside parallel_reduce(...)
#     bin_count_file = pll_reduce.parallel_reduce(
#         ngs_mcount_reduce, temp_count_filenames, processes= 16)
#     
#     #reduce should have made this list have a length one 
#     # (one final bin count file)
#     assert len(bin_count_file) == 1
#     
#     #name for the bin count file, also add it to a list of file names
#     bin_count_files[bin] = (
#         output_prefix +'.%s' % (bin) + '.counts')
#     
#     #copy the bin_count_file out of the temp directory and get rid of
#     #the original
#     shutil.copy(bin_count_file[0],bin_count_files[bin])
#     
#     return bin_count_files

# ============================================================================
# = MULTIPLE THREAD READ COUNT MAP/REDUCE
# ============================================================================

# def ngs_mcount_reduce(count_file1, count_file2):
#     ''' take two count files and merge them into one by combining counts'''
#     
#     count_file1 = iter(open(count_file1, 'rU'))
#     count_file2 = iter(open(count_file2, 'rU'))
#     
#     new_count_file = tempfile.NamedTemporaryFile(delete= False)
#     
#     #print 'Reduce file: ',new_count_file.name
#     
#     line1 = count_file1.next()
#     line2 = count_file2.next()
#     
#     while True:
#         try:            
#             (count1, seq1) = line1.split()
#             (count2, seq2) = line2.split()
#             cmpval = cmp(seq1, seq2)
#             
#             
#             if cmpval == 0:
#                 merged_count = int(count1) + int(count2)
#                 print >> new_count_file, ' '.join((str(merged_count), seq1))
#                 line1 = count_file1.next()
#                 line2 = count_file2.next()
#                 
#             elif cmpval < 0 :
#                 print >> new_count_file, line1.lstrip().rstrip()
#                 line1 = count_file1.next()
#                 
#             elif cmpval > 0:
#                 print >> new_count_file, line2.lstrip().rstrip()
#                 line2 = count_file2.next()
#             
#             continue
#                 
#         except StopIteration:
#             #either both files ended, or 
#             #one of the files has ended before the other,
#             # write the remainder of both files to the merged one
#             
#             #print 'Reached stop iteration: ',new_count_file.name
#             
#             remaining1 = map(lambda line: line.lstrip().rstrip(), 
#                 list(count_file1))
#             remaining2 = map(lambda line: line.lstrip().rstrip(), 
#                 list(count_file2))
#                 
#             if cmpval > 0:
#                 print >> new_count_file, '\n'.join(
#                     [line1.lstrip().rstrip()] + remaining1)
#                     
#             elif cmpval < 0:
#                 print >> new_count_file, '\n'.join(
#                     [line2.lstrip().rstrip()] + remaining2)
#                     
#             elif cmpval == 0:
#                 if remaining1: print >> new_count_file, '\n'.join(remaining1)
#                 if remaining2: print >> new_count_file, '\n'.join(remaining2)
#                 
#             break
#         
#     count_file1.close()
#     count_file2.close()
#     new_count_file.close()
# 
#     return new_count_file.name
# 
# def ngs_mcount_map(seq_file, **kwargs):
#     ''' for multi-thread counting, use the shell to count unique reads per
#         tile/bin, then return a file-like obj with the counts for that bin.'''
#     
#     seq_file = util.open_zcat(seq_file)
#     count_file = tempfile.NamedTemporaryFile(delete=False)
#     
#     #print 'Count File to map to: ', count_file.name
#     
#     command = 'grep -Pi \'^[ATGC]+$\' | sort | uniq -c '
#     
#     subprocess.call(command,
#         stdout=count_file,
#         stdin=seq_file,
#         shell=True,
#         executable='/bin/bash') 
#     
#     return count_file.name
# 
#     
# def ngs_count_reads(merged_fq, output_prefix, **kwargs):
#     '''given aligned SeqPrep metareads, 
#         make hash of all sequences in order to count occurences of each
#         return a list of SeqRecord objects with count annotations'''
#         
#     #counter.py taken from:
#     #http://docs.python.org/library/collections.html#collections.Counter
#     
#     print >> sys.stderr, 'Creating count of unique sequences...'
#     sys.stderr.flush()
#     
#     if kwargs['no_manual_trim']:
#         fq_merged_seqs = SeqIO.parse(
#             gzip.open(merged_fq, "rb"), "fastq-sanger")
#     else:
#         fq_merged_seqs = []
#         trimmed_seq_file = gzip.open(merged_fq, 'rb')
#         for line in trimmed_seq_file:
#             fq_merged_seqs.append(SeqRecord(
#                 seq=Seq(line.rstrip(),IUPAC.unambiguous_dna)))
#     
#     if kwargs['revcom']:
#         fq_counts = counter.Counter(
#             [record.seq.reverse_complement().tostring()
#                 for record in fq_merged_seqs])
#         
#     else:
#         fq_counts = counter.Counter(
#             [record.seq.tostring() for record in fq_merged_seqs])
#     
#     return fq_counts