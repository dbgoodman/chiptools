import tempfile
import subprocess
import sys
import os
import re
import tempfile
import glob
import time
import pdb
import itertools

from multiprocessing import Pool
from collections import defaultdict
from StringIO import StringIO
from operator import truediv, itemgetter

import util
import chunk_fasta
import config

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

def call_seqprep(output_prefix, fq_read1_fn, fq_read2_fn,
    adapter1, adapter2, **kwargs):

    #generate file names
    files = {
            'fq1_out_fn': output_prefix+'.R1.sp.fq.gz',
            'fq2_out_fn': output_prefix+'.R2.sp.fq.gz',
            'fq1_disc_fn': output_prefix+'.R1.disc.fq.gz',
            'fq2_disc_fn': output_prefix+'.R2.disc.fq.gz',
            'merged_out_fn': output_prefix+'.M.fq.gz'}

    #skip if the option is turned on.
    if kwargs['skip_finished'] and os.path.exists(files['merged_out_fn']):
        return files

    #make command
    command = '''
        %(seqprep_path)s
        -f %(fq_read1_fn)s
        -r %(fq_read2_fn)s
        -1 %(fq1_out_fn)s
        -2 %(fq2_out_fn)s
        -3 %(fq1_disc_fn)s
        -4 %(fq2_disc_fn)s
        -s %(merged_out_fn)s
        -A %(adapter1)s -B %(adapter2)s
        -X 1 -g -L 5''' % (dict(locals().items() + files.items() + config.bin_paths.items()))

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
    stats['command'] = command
    files['stats'] = stats

    return files

def sp_stats(fq_tree, output_prefix, **kwargs):
    ''' statistics file for each seqprep tile file, pairs merged, trimmed,
    discarded, CPU times, etc.'''

    bin_totals = defaultdict(int)

    if kwargs['skip_finished'] and os.path.exists(output_prefix+'.spstats'):
        with open(output_prefix+'.spstats') as spstats_file:
            #discard header:
            spstats_file.next()

            for line in spstats_file:
                vals = line.split()
                bin_totals[int(vals[0])] += int(vals[3])
        return bin_totals
    else:
        with open(output_prefix+'.spstats', 'w') as spstats_file:

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
                for tile, tdict in sorted(fq_tree[bin].items()):
                    ppro = tdict['stats']['Pairs Processed']
                    pmrg = tdict['stats']['Pairs Merged']
                    pada = tdict['stats']['Pairs With Adapters']
                    pdsc = tdict['stats']['Pairs Discarded']
                    cput = tdict['stats']['CPU']

                    bin_totals[int(bin)] += pmrg

                    #print line:
                    print >> spstats_file, '\t'.join(map(str,[
                        bin,tile,ppro,
                        pmrg, pmrg/float(ppro),
                        pada, pada/float(ppro),
                        pdsc, pdsc/float(ppro),
                        cput]))
    return bin_totals

def ngs_make_fqtree(fq_reads, regex_str, **kwargs):

    '''load/extract fastq/fasta files:
        we just want to check that the fastq files look right.
        We also want to group the reads by tile, read number, and bin into
        a dict tree that goes bin - tile - read'''

    fq_tree = defaultdict(lambda: defaultdict(dict))

    #create a tree of all the fastq files
    for fq_file in fq_reads:
        fq_components = re.match(regex_str, fq_file.name)

        if not fq_components:
            raise IOError("%s does not match supplied regex %s" % (fq_file.name, regex_str))

        fq_read = fq_components.group('read')
        fq_bin = fq_components.group('bin')
        fq_tile = fq_components.group('tile')

        fq_tree[fq_bin][fq_tile][fq_read] = fq_file

    #check the tree. print out the number of bins, the number of tiles, and
    #make sure that each bin has the same tiles and each tile
    #has only two reads, 1 and 2.

    #REMOVED: check that tiles are the same
    #(NOT DOING NOW, THIS FOR HISEQ)
    #all_tiles = [frozenset(bin.keys()) for bin in fq_tree.values()]
    #if not len(set(all_tiles)) <= 1:
    #    raise( ValueError('Tiles are not the same among bins: '
    #        + print_fq_tree(fq_tree)))

    #check that each tile has only two reads, called 1 and 2
    all_reads = [
        frozenset([frozenset(bin[tile].keys()) for tile in bin.keys()])
            for bin in fq_tree.values()]
    if not len(set(all_reads)) <= 1:
        raise( ValueError('Not all tiles have the same two read files: '
            + print_fq_tree(fq_tree)))

    if not set(all_reads).pop() == frozenset([frozenset(['1', '2'])]):
        raise( ValueError('All reads must be named 1 and 2: '
            + print_fq_tree(fq_tree)))

    return fq_tree

def ngs_maptiles(output_prefix, fq_tree, **kwargs):
    '''
    Send trim and sort jobs (seqprep) and merge and count jobs to workers
    in a pool. Create temporary files per bin of unique counted strings.
    '''

    print >> sys.stderr, 'Trimming, sorting, and counting tiles...\n'

    pool = Pool(None)
    ts_jobs = defaultdict(dict) #trimsort jobs
    mc_jobs = {} #mergecount jobs

    bin_files = {}

    for bin in fq_tree.keys():
        for tile, files in fq_tree[bin].items():
            fq_read1_fn = files['1'].name
            fq_read2_fn = files['2'].name

            #tile and bin prefix
            t_b_prefix = output_prefix+'.%s.%s' % (bin, tile)

            ts_jobs[bin][tile] = pool.apply_async(ngs_trimsort,
                (t_b_prefix, fq_read1_fn, fq_read2_fn,),
                kwargs, callback= files.update)
            # ngs_trimsort(t_b_prefix, fq_read1_fn, fq_read2_fn,**kwargs)

    while True:
        t_done = 0
        t_tot = 0

        #if all tiles for a bin have completed, merge and count them
        for bin, tiles in ts_jobs.items():

            t_done += sum([ts.ready() for ts in tiles.values()])
            t_tot += len(tiles)

            #if the job for this bin has been started, then skip
            if bin in mc_jobs: continue

            #if all tiles are complete:
            if all([job.ready() for job in tiles.values()]):

                try: assert job.successful()
                except AssertionError: job.get()

                file_sets = fq_tree[bin].values()
                mc_jobs[bin] = pool.apply_async(ngs_mergecount,
                    (file_sets, bin), kwargs,
                    callback= bin_files.update)

        b_done = sum([mc.ready() for mc in mc_jobs.values()])
        b_tot = len(fq_tree.keys())

        print >> sys.stderr, '\t%2d / %2d tile jobs completed.' % (
                t_done, t_tot),
        print >> sys.stderr, '\t%2d / %2d bin jobs completed.\r' % (
                b_done, b_tot),
        sys.stderr.flush()

        if t_done == t_tot and b_done == b_tot:
            for job in mc_jobs.values():
                try: assert job.successful()
                except AssertionError: job.get()
            print >> sys.stderr, '\n'
            break

        else:
            time.sleep(1)

    return bin_files, fq_tree

def ngs_trimsort(output_prefix, fq_read1_fn, fq_read2_fn,
    **kwargs):
    '''Steps to perform with bash pipes:
    1. call seqprep and get a trimmed merged FQ file
    2. convert FQ to seq (FA w/ no headers), manually trim reads with adapters
    3. sort, gzip, and return files in files{} dict
    '''

    files = call_seqprep(output_prefix, fq_read1_fn, fq_read2_fn, **kwargs)

    files['trimmed_out_fn'] = output_prefix+'.MT.seq.gz'

    #skip if the option is turned on.
    if kwargs['skip_finished'] and os.path.exists(
        files['trimmed_out_fn']):
        return files

    #either just grab sequences and sort, or do manual trim:

    if kwargs['no_manual_trim']:
        #do not look for adapters, just remove FASTA headers
        trim_cmd = '%(zcat_path)s %(merged_out_fn)s | %(grep_path)s -P \'^[NATGC]+$\'' % (
            dict(kwargs.items() + files.items() + config.bin_paths.items()))
    else:
        #cull remaining sequences whose adapters were not trimmed
        trim_cmd = ' '.join([
            'cat <(%(zcat_path)s %(merged_out_fn)s | %(grep_path)s -P \'^[NATGC]+$\'',
                 '| %(grep_path)s -Pv \'(?<=^|%(rs2)s)[NATGC]+?(?=$|%(rs1)s)\')',
                 '<(%(zcat_path)s %(merged_out_fn)s | %(grep_path)s -P \'^[NATGC]+$\'',
                 '| %(grep_path)s -Po \'(?<=^|%(rs2)s)[NATGC]+?(?=$|%(rs1)s)\')']) % (
                    dict(kwargs.items() + files.items() + config.bin_paths.items()))

    #we might be taking the reverse complement of all sequences first...
    if kwargs['revcom']:
        trim_cmd += ' | perl -ne \'chomp; $_ =~ tr/ACGTacgt/TGCAtgca/;'
        trim_cmd += ' print reverse($_)."\\n";\''

    # if kwargs['RNA']:
    #     trim_cmd += ' | perl -ne \'length($_) > 3 && print substr $_, 2;\''

    #filter on minimum read length if specified at command line
    if 'minlen' in kwargs:
        trim_cmd += ' | perl -ne \'(length($_) >= %(minlen)s) && print;\'' % (
            dict(kwargs.items()))

    trim_cmd += ' | sort | gzip > %(trimmed_out_fn)s' % (
        dict(kwargs.items() + files.items()))

    subprocess.call(trim_cmd, shell=True, executable='/bin/bash')

    assert len(files) > 0
    return files

def ngs_mergecount(file_sets, bin, **kwargs):
    '''Steps to perform with bash pipes:
    1. merge and count a list of trimsort results with sort -m and uniq -c
    3. gzip and return files in files{} dict
    '''

    seq_files = util.open_zcat(
        [files['trimmed_out_fn'] for files in file_sets])

    count_file = tempfile.NamedTemporaryFile(delete=False)

    command = ' '.join(['sort -m | uniq -c | perl -ne',
        '\'split && print "$_[1] %(bin)s $_[0]\n"\' | gzip']) % locals()

    subprocess.call(command,
        stdout=count_file,
        stdin=seq_files,
        shell=True,
        executable='/bin/bash')

    return {bin:count_file.name,}

def final_counts(output_prefix, bin_files, **kwargs):
    '''now we want to open file objects for all bins and create a fasta file
    with all the counts in a tabbed list in sorted bin order. each record
    will look like:
    >SEQ_NUM (\t) tot_val (\t) bin1_val (\t) bin2_val (\t) ... (\t) bin12val
    SEQUENCE
    '''

    print >> sys.stderr, 'Merging bin counts...\n'

    all_bins_sorted = subprocess.Popen('sort',
        stdout=subprocess.PIPE,
        stdin=util.open_zcat(bin_files.values()),
        shell=True,
        executable='/bin/bash')

    #open the fa file
    counts_file = open(output_prefix+'.counts.fa', 'w')
    seq_num = 0

    #bins must be consecutive!
    counts = [0] * len(bin_files)

    #get first line
    seq, bin, count = all_bins_sorted.stdout.next().split()
    counts[int(bin)-1] += int(count)
    prev_seq = seq

    line = ">%(seq_num)d\t%(tot_count)d\t%(bin_str)s\n%(prev_seq)s"

    while 1:
        try:
            #get next line
            seq, bin, count = all_bins_sorted.stdout.next().split()

            #update counts
            if seq == prev_seq:
                 counts[int(bin)-1] += int(count)

            #print line
            else:

                tot_count = sum(counts)
                bin_str = '\t'.join(map(str, counts))

                if tot_count >= kwargs['tot_count']:
                    print >> counts_file, line % locals()
                    seq_num += 1

                prev_seq = seq

                #populate the counts array with this first value
                counts = [0] * len(bin_files)
                counts[int(bin)-1] += int(count)

        except StopIteration:
            break

def ngs_cluster(fq_files, regex_str, output_prefix, **kwargs):

    fq_tree = ngs_make_fqtree(fq_files, regex_str, **kwargs)

    bin_files, fq_tree = ngs_maptiles(output_prefix, fq_tree, **kwargs)

    bin_totals = sp_stats(fq_tree, output_prefix, **kwargs)

    final_counts(output_prefix, bin_files, **kwargs)
