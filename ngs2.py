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
        /opt/SeqPrep.dbg/SeqPrep
        -f %(fq_read1_fn)s
        -r %(fq_read2_fn)s
        -1 %(fq1_out_fn)s
        -2 %(fq2_out_fn)s
        -3 %(fq1_disc_fn)s
        -4 %(fq2_disc_fn)s
        -s %(merged_out_fn)s
        -A %(adapter1)s -B %(adapter2)s
        -X 1 -g -L 5''' % (dict(locals().items() + files.items()))

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
    discarded, CPU times, etc. Also retur'''

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
        trim_cmd = 'zcat %(merged_out_fn)s | grep -P \'^[NATGC]+$\'' % (
            dict(kwargs.items() + files.items()))
    else:
        #cull remaining sequences whose adapters were not trimmed
        trim_cmd = ' '.join([
            'cat <(zcat %(merged_out_fn)s | grep -P \'^[NATGC]+$\'',
                 '| grep -Pv \'(?<=^|%(rs2)s)[NATGC]+?(?=$|%(rs1)s)\')',
                 '<(zcat %(merged_out_fn)s | grep -P \'^[NATGC]+$\'',
                 '| grep -Po \'(?<=^|%(rs2)s)[NATGC]+?(?=$|%(rs1)s)\')']) % (
                    dict(kwargs.items() + files.items()))

    #we might be taking the reverse complement of all sequences first...
    if kwargs['revcom']:
        trim_cmd += ' | perl -ne \'chomp; $_ =~ tr/ACGTacgt/TGCAtgca/;'
        trim_cmd += ' print reverse($_)."\\n";\''

    if kwargs['RNA']:
        trim_cmd += ' | perl -ne \'length($_) > 3 && print substr $_, 2;\''

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

def ngs_usearch_wdb(output_prefix, library_fn, **kwargs):
    usearch_wdb_cmd = ' '.join(['usearch --makewdb %(library_fn)s',
        '--output %(output_prefix)s.wdb -w 4 --slots 400000003']) % locals()

    subprocess.call(usearch_wdb_cmd, stdout= subprocess.PIPE,
        stderr=subprocess.PIPE, shell=True)

def ngs_usearch(output_prefix, input_fn, bin_totals, **kwargs):

    #global values for parse_line(...)
    #-------------------------------------------------------------------------
    #pre-calculate bin and header associated strings and arrays
    nbins = len(bin_totals)
    if kwargs['RNA']:
        num_orig_fields = (6 + nbins)
    else:
        num_orig_fields = (6 + nbins)

    #for now, specify the bin means by hand:
    bin_means = [35.35533906, 1591.775738, 2581.230133, 4186.047897,
        6788.703484, 11008.82623, 17852.78592, 28952.36519, 46952.04117,
        76141.41273, 123477.743, 200242.2658][:nbins]

    #precompile regular expression matches for the cigar string
    ins_pat = re.compile('(\d*)I')
    mat_pat = re.compile('(\d*)M')
    sum_ins = lambda cigar: sum(map(int,
        [(i or 1) for i in ins_pat.findall(cigar)])) or 0
    sum_mat = lambda cigar: sum(map(int,
        [(i or 1) for i in mat_pat.findall(cigar)])) or 0

    #-------------------------------------------------------------------------
    def parse_line(line):
        '''
        This function will be called one per printed line. It takes the values
        generated by USEARCH and from them calculates additional values more
        useful to me, including additional match-based statistics and read
        expression values.
        ---
        initial fields: (fields that USEARCH generates)
            0:query, n+2:target ,n+3:id, n+4:caln, n+5:pv, n+6:pvz, n+7:gaps
                pv: matching columns #qlen if RNA
                pvz: nonmatching columns
                gaps: ins + del

                if RNA: (not printing these at the moment)
                pv is qlen
                tloz: zero-based target start
                thiz: zero-based target end
                qloz: zero-based query start
                qhiz: zero-based query end

            final output fields:
            old: (just copy these from USEARCH)
                0:query
                1:total_reads
                2 .. n+2:bins
                n+2:target
                n+3:id
                n+4:caln
                n+5:nmatch
                n+6:mismatch
                (if rna: n+7 through n+10: tloz, thiz, qloz, qhiz)

            new: (calculate these now and append them in order)
                ins, del, lev, adj_total
            adj_total: scaled and adjusted number of reads by mean bin value

        '''

        field_vals = line.split()
        new_vals = field_vals[0:num_orig_fields]
        perfect = (new_vals[0].lstrip('>') in kwargs['perfects'] and
            len(kwargs['perfects'][new_vals[0].lstrip('>')]) > 0)

        num_perfect = 0

        #if this read has already been perfectly matched without UCLUST, then
        #generate its old vals automatically
        if perfect:
            #remove extra values from UCLUST if they exist
            if kwargs['RNA']:
                #seq_dict contains lib_id, and unmatched prefix, and suffixes
                seq_dict = kwargs['perfects'][new_vals[0].lstrip('>')][0]

                #get lengths of lib_seq, prefix, suffix
                seq_id = seq_dict['lib_id']
                seq_len = seq_dict['seq_len']
                num_perfect = len(kwargs['perfects'][new_vals[0].lstrip('>')])
                lib_len = kwargs['lib_len'][seq_id]
                seq_start = len(seq_dict['start'])
                after_ins = lib_len - seq_start - seq_len

            else:
                seq_id = kwargs['perfects'][new_vals[0].lstrip('>')]
                seq_len = kwargs['lib_len'][seq_id]
            new_vals = new_vals[0:(nbins+2)]
            #target---
            new_vals.append(seq_id)
            #id---
            new_vals.append('100.0')
            #cigar---
            if kwargs['RNA']:
                cigar = []
                if seq_start > 0: cigar.append(str(seq_start)+'I')
                cigar.append(str(seq_len)+'M')
                if after_ins > 0: cigar.append(str(after_ins)+'I')
                new_vals.append(''.join(cigar))
            else:
                new_vals.append(str(seq_len)+'M')
            #match---
            new_vals.append(str(seq_len))

            #Remaining fields:---
            #RNA: mm = 0, ins = before and after, del = 0, lev = 0
            if kwargs['RNA']: new_vals.extend([0, after_ins+seq_start, 0,0])
            #DNA: mm, ins, del, lev are all 0
            else: new_vals.extend([0]*4)

            # #if RNA, then query and target start and end are also 0
            # if kwargs['RNA']:
            #     new_vals.extend([0]*4)

        #if this read does not match, then blank out its match-related fields
        elif len(field_vals) < num_orig_fields:
            new_vals[0] = new_vals[0].lstrip('>')
            #fill in match associated values with '*'
            new_vals.extend(['NA'] * (uclust_raw_num_cols+1))

        #calculate for reads that match to a designed sequence
        else:
            #mismatch: pv - pvz
            #new_vals.append(str(sum_mat(field_vals[nbins+4]))+'-'+new_vals[-1])
            new_vals.append(sum_mat(field_vals[nbins+4]) - int(new_vals[-1]))
            #ins: use regex to parse insertions
            new_vals.append(sum_ins(field_vals[nbins+4]))
            #del: subtract ins from gaps (ins = dels + gaps)
            new_vals.append(int(field_vals[nbins+7]) - int(new_vals[-1]))
            #lev: mismatch+del+ins
            new_vals.append(sum(map(int, new_vals[-3:])))

        #adj_total: (for both matching and non-matching reads)
        #==========================================================

        #divide the number of reads by the total reads for that bin
        adj_totals = [truediv(float(a), b)
            for a, b in zip(field_vals[2:nbins+2], bin_list)]

        #sum the adjusted read totals
        adj_sum = sum(adj_totals)

        #divide the reads in each bin by the adjusted total reads to
        #get an adjusted read pct per bin
        adj_totals = [i/adj_sum for i in adj_totals]

        #multiply each adjusted pct per bin by the mean expr of the bin
        adj_val = sum([i*j for i,j in zip(adj_totals, bin_means)])
        new_vals.append(adj_val)

        if not kwargs['RNA']:
            if len(field_vals) < num_orig_fields or perfect:
                new_vals.append('NA')
            else:
                new_vals.append(field_vals[-1]) #alts
        elif kwargs['RNA']:
            if perfect and num_perfect > 1:
                new_vals.append(
                    '; '.join(map(itemgetter('lib_id'),
                        kwargs['perfects'][new_vals[0].lstrip('>')])))
            else:
                new_vals.append('NA')


        # if new_vals[0] == '6591':
        #     pdb.set_trace()

        return '\t'.join(map(str,new_vals))+'\n'
    #-------------------------------------------------------------------------

    #turn the bin totals dict into a list
    bin_list = []
    [bin_list.append(total) for bin, total in sorted(bin_totals.items())]

    usearch_out_file = tempfile.NamedTemporaryFile(delete=False)
    usearch_out_fn = usearch_out_file.name
    usearch_out_file.close()

    if kwargs['RNA']:
        #FOR RNA:
        usearch_cmd = ' '.join([
            'usearch --query %(input_fn)s -w 8',
            '--wdb %(output_prefix)s.wdb --maxaccepts 0 --maxrejects 0',
            '--id 1 --weak_id 0.5 --minlen 25 --local',
            '--idsuffix 20',
            #'--gapopen 0TE/0QE/*QI/*TI',
            '--rightjust',
            '--mismatch -4',
            #'--nofastalign',
            '--iddef 2',
            '--userout %(usearch_out_fn)s',
            '--userfields',
            'query+target+id+caln+ql+pvz+gaps']) % locals()
        uclust_raw_num_cols = 7

        #TESTING:
        #switch from: 'query+target+id+caln+pv+pvz+gaps+tlowz+thiz+qloz+qhiz']) % locals()
        #switch to: 'query+target+id+qrow+pv+pvz+gaps+trow+thiz+qloz+qhiz'
        #cut -f5,6,7,11 \
        #    scratch/dbg/ecre/ct/203_hsrna_small/203_hsrna_small.ustats \
        #    | perl -ne 's/(^[\w-]+\s+[\d.]+)\s+([-\w]+)\s+([-\w]+)/>$1\nQ:
        #        $2\nT: $3/ && print' \
        #    > /scratch/dbg/ecre/ct/203_hsrna_small/203_hsrna_small.falign

    else:
        #FOR DNA (default):
        usearch_cmd = ' '.join([
            'usearch --query %(input_fn)s -w 6',
            '--wdb %(output_prefix)s.wdb --maxaccepts 1 --maxrejects 50',
            '--id 1 --weak_id 0.5 --global --minlen 25 --iddef 1 -w 4',
            '--userout %(usearch_out_fn)s',
            '--userfields query+target+id+caln+pv+pvz+gaps']) % locals()
        uclust_num_cols = 7

    #print usearch_cmd

    subprocess.call(usearch_cmd, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, shell=True)

    ustat_file = tempfile.NamedTemporaryFile(delete=False)
    input_file = open(input_fn, 'r')

    id_lines = []
    curr_id = -1
    best_pct = 0
    best_lines = []
    fa_lines = []

    #-parsing of USEARCH OUTPUT-----------------------------------------------
    # due to a bug, USEARCH prints out every possible hit for nrejects or
    # naccepts. I am going to pull out the best hit among all of those found.
    # I am also adding in unmatched reads in order of sequence as well.

    for line in open(usearch_out_fn,'r'):
        line_vals = line.split()

        #skip if header (will have 7 lines)
        if len(line_vals) == uclust_raw_num_cols: continue
        line_id, line_pct = [line_vals[i] for i in (0,nbins+3)]

        #-add in unmatched lines----------------------------------------------
        #get the corresponding FASTA line. If it is missing, add it as an
        # unmatched record
        while 1:
            try:
                if line_id == curr_id: break
                fa_line = input_file.next()
                fa_vals = fa_line.split()
                fa_id = fa_vals[0].lstrip('>')

                if fa_id == line_id:
                    input_file.next() #skip the sequence line
                    break

                elif fa_id != line_id:
                    fa_lines.append(fa_line)
                    input_file.next() #skip the sequence line
                    continue
            except StopIteration:
                break
        #---------------------------------------------------------------------
        line_pct = float(line_pct)
        if curr_id == line_id:
            if line_pct > best_pct:
                best_pct = line_pct
                best_lines = [line]
            elif line_pct == best_pct:
                best_lines.append(line)
            else:
                continue
        else:
            if best_lines:
                #if there are multiple designed sequences that match this
                #read to the same level of identity, then store them in the
                #'alts' column separated by ';'s
                if len(best_lines) > 1:
                    best_hits = [bl.split()[nbins+2] for bl in best_lines[1:]]
                    best_hits = ('\t'+';'.join(best_hits))
                else:
                    best_hits = '\tNA'

                out_lines = ''.join(
                    map(parse_line, [best_lines[0]+best_hits] + fa_lines))

                #print out_lines
                ustat_file.write(out_lines)

                best_lines = [line]
                fa_lines = []

            curr_id = line_id
            best_pct = line_pct
        #--end of USEARCH parsing---------------------------------------------

    os.remove(usearch_out_fn)
    os.remove(input_fn)
    return ustat_file.name

def RNA_find_perfect(lib_str,rec_num,rec_seq):
    '''simple paralllelized function that uses a quick RE to find
    all instances of the record string as a substring of a longer
    library member'''

    perfect_matches = [m.groupdict() for m in re.finditer(
        (r'\|(?P<lib_id>[\w\-*]+)\|(?P<start>\w*)' + rec_seq), lib_str)]
    [pm.update(seq_len= len(rec_seq)) for pm in perfect_matches]
    return {rec_num : perfect_matches}


def ngs_parallel_usearch(output_prefix, bin_totals, **kwargs):
    ''' I am taking my unique output fasta file and searching for best
        hits with grep followed by USEARCH. Grep will find perfect matches
        only, which USEARCH might miss due to its heuristics.
    '''


    print >> sys.stderr, 'Searching for perfect matches...\n'

    #first we have to make a string dict from the library file:
    lib_fh = open(kwargs['lib_records'].name)
    rec_fh = open(output_prefix + '.counts.fa')


    lib_dict = dict()
    kwargs['lib_len'] = dict()
    perfect_matches = dict()
    lib_str = ''

    try:
        while 1:

            lib_id = lib_fh.next().lstrip('>').rstrip()
            lib_seq = lib_fh.next().rstrip().upper()
            lib_dict[lib_seq] = lib_id
            lib_str += '|' + lib_id + '|' + lib_seq
            kwargs['lib_len'][lib_id] = len(lib_seq)
    except StopIteration:
        pass

    if kwargs['RNA']:
        rna_jobs = []
        pool = Pool()

    rec_tot = 0

    for line in rec_fh:
        rec_tot += 1
        rec_num = line.lstrip('>').split('\t',1)[0]
        rec_seq = rec_fh.next().rstrip().upper()

        if kwargs['RNA']:
            #RNA perfect match requires parallel loop-------------------------

            rna_jobs.append(pool.apply_async(RNA_find_perfect,
                 (lib_str, rec_num, rec_seq),
                 callback= perfect_matches.update))
            # perfect_matches[rec_num] = RNA_find_perfect(
            #     lib_str,rec_num, rec_seq)

        else:
            #DNA MATCH is just a simple search--------------------------------
            found = lib_dict.get(rec_seq, False)
            if found:
                perfect_matches[rec_num] = found

    if kwargs['RNA']:
        #wait for all parallel RNA jobs to finish-------------------------
        while True:
            rna_done = sum([rna.ready() for rna in rna_jobs])

            update_str = '\t%d / %d RNA match jobs completed.\r'
            print >> sys.stderr, update_str % (rna_done, rec_tot),
            sys.stderr.flush()

            #if all tiles are complete:
            if rna_done == rec_tot:
                for job in rna_jobs:
                    try: assert job.successful()
                    except AssertionError: job.get()
                break
            else: time.sleep(1)

        pool.close()
        pool.join()

    missing_fh = open(output_prefix + '.missing.txt','w')

    if kwargs['RNA']:
        missing_perfects = set(lib_dict.values()) - set(
            itertools.imap(itemgetter('lib_id'), itertools.chain(
                *perfect_matches.values())))
    else:
        missing_perfects = set(lib_dict.values()) - set(
            perfect_matches.values())

    #this will be passed to the UCLUST functions so that they are skipped
    #and written to the ustats file differently
    kwargs['perfects'] = perfect_matches


    for lib_name in missing_perfects:
        print >> missing_fh, lib_name

    print >> sys.stderr, (
        '\t %d missing perfect library members.\n' % (len(missing_perfects)))
    print >> sys.stderr, 'Splitting FASTA unique file into chunks...\n'

    fa_chunk_fns = chunk_fasta.chunk_file(
        output_prefix + '.counts.fa', split_type='pieces', mp=False)
    print >> sys.stderr, 'Clustering with USEARCH...\n'

    #start multiprocessing pool
    pool = Pool()
    output_files = []
    us_jobs = []

    for fa_chunk_fn in fa_chunk_fns:
        # us_jobs.append(pool.apply_async(ngs_usearch,
        #      (output_prefix, fa_chunk_fn, bin_totals), kwargs,
        #      callback= output_files.append))
        output_files.append(
             ngs_usearch(output_prefix, fa_chunk_fn, bin_totals, **kwargs))


    while True:
        us_done = sum([us.ready() for us in us_jobs])
        us_tot = len(fa_chunk_fns)

        print >> sys.stderr, '\t%d / %d usearch jobs completed.\r' % (
                 us_done, us_tot),
        sys.stderr.flush()

        #if all tiles are complete:
        if us_done == us_tot:
            for job in us_jobs:
                try: assert job.successful()
                except AssertionError: job.get()
            break

        else:
            time.sleep(1)

    pool.close()
    pool.join()

    ustats_file = open(output_prefix+'.ustats', 'w')
    bin_names = ['bin.'+str(i) for i in range(1,len(bin_totals)+1)]
    field_names = (['uniq_id', 'counts'] + bin_names
        + ['target','id','caln','match','mismatch',
           'ins','del','lev','score','alts'])
    ustats_file.write('\t'.join(field_names)+'\n')
    ustats_file.flush()
    subprocess.call('cat '+' '.join(output_files), #+' | sort -nrk2',
        stdout= ustats_file, shell=True, executable='/bin/bash')

    #finally, delete the temporary files
    map(os.remove, output_files)
    print >> sys.stderr, '\n'

def ngs_cluster(fq_files, regex_str, output_prefix, **kwargs):

    fq_tree = ngs_make_fqtree(fq_files, regex_str, **kwargs)

    bin_files, fq_tree = ngs_maptiles(output_prefix, fq_tree, **kwargs)

    bin_totals = sp_stats(fq_tree, output_prefix, **kwargs)

    final_counts(output_prefix, bin_files, **kwargs)

    #skip uclust steps if the 'no_uclust' flag is set
    if kwargs['no_uclust']:
        return

    ngs_usearch_wdb(output_prefix, kwargs['lib_records'].name, **kwargs)

    ngs_parallel_usearch(output_prefix, bin_totals, **kwargs)
