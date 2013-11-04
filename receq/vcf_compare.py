import vcf
import re
import glob
import sys
import subprocess
import pdb
import os
import time
import itertools
import pickle

import numpy as np

from operator import itemgetter
from multiprocessing import Pool
from shove import Shove
from StringIO import StringIO
from collections import defaultdict, namedtuple

from pandas import *

#Reader class can't be pickled without these lines:

vcf.parser.Info = vcf.parser._Info
vcf.parser.Format = vcf.parser._Format
vcf.parser.Filter = vcf.parser._Filter
vcf.parser.Call = vcf.parser._Call
vcf.model.CallData = vcf.model._CallData
#
#Inputs:
# Reference Genome (FASTA)
# BAM File (.bam)

#Steps:
# Call haploid and diploid with Freebayes (freebayes.sh)
# Call diploid with GATK (gatk.sh)
# Call effects with snpEff (snpeff.sh)
# Get records for each using a pipe
# Combine records into one
# Write it out to a new VCF using a writer
# Generate proxy statistics: "What information do I want..."
# Also write a more 'readable' tab delimited file

#What information do I want:
# * does this SNP have a potential phenotype (in a gene, non-synonymous, etc)
#    run snpEFF on every position (pipe a writer to snpeff.sh) 
#    and get the EFF field, specifically the type, severity, and gene
#    EFF_TYPE(SEVERITY|...|...|...|GENE)
#    eg: EFF=CODON_DELETION(MODERATE||gcgaaagcggtg/gcg|AKAV414A|deoA|...)
# * is it in a duplicated region
#    a proxy for this is mean mapping quality of ref vs alt alleles:
#     MQMR (reference), MQM (alternate)
# * what is coverage, and what percentage of reads support it
# * is the fwd/reverse read ratio the same for mut vs wt
# * which callers find it? i.e. freebayes hap, freebayes dip, gatk, breseq
# * print out a tview-style alignment

# vcf_filenames = sys.argv[1:]
# #vcf_fn_regex = re.compile(r'.*/([\.\w]+)\.vcf')
# vcf_fn_regex = re.compile(r'.*/([\w]+).([\.\w]+)\.vcf')
# 
# vcf_readers = {}
# vcf_wstreams = {}
# 
# vcf_records = Read
# 
# for vcf_fn in vcf_filenames:
#         
#     vcf_fn_match = vcf_fn_regex.match(vcf_fn)
#     sample, descriptors = vcf_fn_match.groups()
#     descriptors = descriptors.split('.')
#     
#     vcf_reader = vcf.Reader(open(vcf_fn, 'rb'))
#     temp_vcf_file = tempfile.NamedTemporaryFile(delete= False)
#     vcf_writer = vcf.Writer(vcf_reader)


def call_snpeff(vcf_process, vcf_outfile='/dev/null', html_outfile=None, 
    file=False, **kwargs):
    
    '''
    call snpeff through the shell on a given vcf input file. optionally
    save the snpeff file (will be redundant with most of the input file)
    and an html output file. Return a vcf reader object.
    '''
    if file:
        vcf_infile = vcf_process
    else:
        vcf_infile = vcf_process.stdout
    
    # if (vcf_outfile and vcf_outfile != '/dev/null' and
    #     os.path.exists(vcf_outfile)):
    #     vcf_process.kill()
    #     return open(vcf_outfile, 'r')
    
    #if we are outputting an html file, then make that line to be added to the
    #shell command below
    if html_outfile:
        html_line = '-s %(html_outfile)s' % (locals())
    else:
        html_line = ''
    
    #create snpeff shell command
    snpeff_cmd = '''
    perl -pe "s/Chromosome/U00096/g" \
    | java -Xmx2G -jar /opt/snpEff/snpEff.jar eff \
        -i vcf -o vcf -c /opt/snpEff/snpEff.config\
        %(html_line)s \
        -ud 50 \
        -no-downstream \
        ecoli |  tee %(vcf_outfile)s''' % dict(
            locals().items() + kwargs.items())
    
    #run snpeff with subprocess
    snpeff_process = subprocess.Popen(
        snpeff_cmd,
        stdin=vcf_infile,
        stdout=subprocess.PIPE,
        stderr=open(os.devnull, 'w'),
        shell=True,
        executable='/bin/bash')
        
    return snpeff_process.stdout

# def get_bedtools_coverage(output_prefix, **kwargs):
# 
#     '''
#     use bedtools to get coverage for a 50 bp window for each genome, and 
#     return the output as a DF
#     '''
#     
#     #if we are outputting an html file, then make that line to be added to the
#     #shell command below
#     if html_outfile:
#         html_line = '-s %(html_outfile)s' % (locals())
#     else:
#         html_line = ''
# 
#     #create snpeff shell command
#     snpeff_cmd = '''
#     perl -pe "s/Chromosome/U00096/g" \
#     | java -Xmx2G -jar /opt/snpEff/snpEff.jar eff \
#         -i vcf -o vcf -c /opt/snpEff/snpEff.config\
#         %(html_line)s \
#         -ud 50 \
#         -no-downstream \
#         ecoli | tee %(vcf_outfile)s''' % dict(
#             locals().items() + kwargs.items())
# 
#     #run snpeff with subprocess
#     snpeff_process = subprocess.Popen(
#         snpeff_cmd,
#         stdin=vcf_infile,
#         stdout=subprocess.PIPE,
#         stderr=open(os.devnull, 'w'),
#         shell=True,
#         executable='/bin/bash')
# 
#     return snpeff_process.stdout

def call_freebayes(bam_filename, fb_params, region= None, **kwargs):
    '''
    call freebayes and return a file stream.
    '''
    
    #print 'FB Region:' + str(region)
    
    
    #if we are debugging or parallelizing, call freebayes only on a small
    #region of the genome
    if region:
        targets_line = "--targets <(echo 'Chromosome %d %d')" % region
    else:
        targets_line = ''
    
    
    bam_filename_line = ' --bam '.join(bam_filename.split())
    
    freebayes_cmd = '''freebayes %(targets_line)s \
        --bam %(bam_filename_line)s \
        --fasta-reference %(fasta_reference_filename)s %(fb_params)s
    ''' % dict(locals().items() + kwargs.items())
    
    #print >> sys.stderr, freebayes_cmd
    
    
    #run snpeff with subprocess
    fb_process = subprocess.Popen(
        freebayes_cmd,
        stdout=subprocess.PIPE,
        stderr=open(os.devnull, 'w'),
        shell=True,
        executable='/bin/bash')
        
    return fb_process

def call_mpileup(bam_filename, region= None, **kwargs):
    '''
    call samtools mpileup and return a file stream.
    '''
    
    #if we are debugging or parallelizing, call freebayes only on a small
    #region of the genome
    
    #print 'MP Region:' + str(region)
    
    if region:
        targets_line = "-l <(echo 'Chromosome %d %d')" % region
    else:
        targets_line = ''
    
    mpileup_cmd = '''samtools mpileup \
        -AuS -C50 -f %(fasta_reference_filename)s \
        %(targets_line)s \
        %(bam_filename)s \
        | bcftools view -vcg -
    ''' % dict(locals().items() + kwargs.items())
    
    #print >> sys.stderr, mpileup_cmd
    
    #run mpileup with subprocess
    mp_process = subprocess.Popen(
        mpileup_cmd,
        stdout=subprocess.PIPE,
        stderr=open(os.devnull, 'w'),
        shell=True,
         executable='/bin/bash')
    
    return mp_process



def write_vcf_parallel(bam_filenames, output_prefix,
    genome= '/scratch/dbg/genomes/mg1655ref.fa', **kwargs):
    
    out_stream = StringIO()
    
    genome_size = 4639675
    chunk_size = 100000
    
    #genome_size = 100000
    #chunk_size = 5000
    
    pool = Pool(12)
    jobs = {} #jobs
    
    chunk_outputs = []
    
    chunks = zip(range(0, (genome_size-chunk_size), chunk_size), 
        range(chunk_size+1,genome_size,chunk_size))
    
    for region in chunks:
                            
        jobs[region] = pool.apply_async(write_vcf,
            (sys.argv[1:], output_prefix, region,), 
            kwargs, callback= chunk_outputs.append)
        # chunk_outputs.append(
        #     write_vcf(bam_filenames, output_prefix, region, **kwargs))
    
    while True:
    
        j_done = sum([j.ready() for j in jobs.values()])
        j_tot = len(chunks)
        
        print >> sys.stderr, '\t%2d / %2d jobs completed.\r' % (
                j_done, j_tot),
        sys.stderr.flush()
    
        if j_done == j_tot:
            for job in jobs.values():
                try: assert job.successful() 
                except AssertionError: job.get()
            print >> sys.stderr, '\n'
            break
            
        else:
            time.sleep(1)
    
    (regions, records_list, metadata) = zip(*chunk_outputs)
    
    ## create a dummy reader class object that holds writer metadata, and 
    ## take the metadata from the first good reader that we find in the chunk
    ## outputs; they should all be the same
    i = 0
    while metadata[i] == None: i += 1
    class DummyReader: pass
    DummyReader.__dict__ = metadata[i]
    
    writer = vcf.Writer(out_stream, DummyReader())
    
    all_records = list(itertools.chain(*records_list))
    
    dict_list = []
    
    for record in all_records:
        writer.write_record(record)
        last_newline = out_stream.getvalue()[:-2].rfind('\n')
        record_line = out_stream.getvalue()[last_newline+1:-2] + '\n'
        
        gts = list(record.REF) + record.ALT
        
        for sample in DummyReader.samples:
            
            sample_call = record.genotype(sample)
            
            if sample_call.data['GT'] == None:
                
                dict_record = {
                    'CHROM': record.CHROM, 
                    'POS': record.POS,
                    'SAMPLE': sample,
                    'REF': record.REF,
                    'TYPE': record.var_subtype,
                    'CALL': None,
                    'VCFLINE': record_line,
                    'RECORD': record
                }
            
            else:    
            
                gt_num = max(map(int,sample_call.gt_nums.split('/')))
            
                dict_record = {
                    'CHROM': record.CHROM, 
                    'POS': record.POS,
                    'SAMPLE': sample,
                    'REF': record.REF,
                    'TYPE': record.var_subtype,
                    'CALL': gt_num,
                    'VAR': sample_call.is_variant,
                    'HET': sample_call.is_het,
                    'BASE': gts[gt_num],
                    'VCFLINE': record_line,
                    'RECORD': record
                }
            
            dict_list.append(dict_record)
    
    df = DataFrame(dict_list,
        columns=['CHROM','POS','SAMPLE','REF','TYPE','CALL','VAR',
                 'HET','BASE','VCFLINE','RECORD'])
    
    pickle.dump(df, open(output_prefix + 'data_frame.p', 'wb'))
    
    out_vcf = open(output_prefix + 'combined_snps.vcf', 'w')
    
    #df = concat(dfs, ignore_index=True)
    
    ## find records where lax freebayes and mp find something, but conserv
    ## freebayes does not
    #bad_positions = df[( (df['fb_recoli_lib3_c31_322A'] == '-') & ((df['fbl_recoli_lib3_c31_322A'] != '-') | (df['mp_recoli_lib3_c31_322A'] != '-') ))]['pos']
    
    ## remove positions where strong freebayes has a snp called but for a 
    ## different record
    #bad_positions = np.unique(bad_positions[np.logical_not(
    #    bad_positions.isin(
    #        df[(df['fb_recoli_lib3_c31_322A'] != '-')]['pos'])) 
    
    #for pos in bad_positions:
    #    print 'Chromosome\t%d' % pos
    
    #playing around:
    
    #concat([df[(df.POS == pos) & (df.VAR == True)][['POS','SAMPLE']] for pos in np.unique(df.POS)]).to_string()
    #df.pivot(index='POS',columns='SAMPLE',values='VAR').to_csv(output_prefix + 'grid.csv')


def write_vcf(bam_filenames, output_prefix, region, 
    genome= '/scratch/dbg/genomes/mg1655ref.fa', **kwargs):
    
    records = []
    
    try:
        reader = vcf.Reader(
            call_snpeff(
                call_freebayes(' '.join(bam_filenames),
                    fb_params=''' \
                        --pvar 0.001 \
                        --ploidy 2 \
                        --min-alternate-fraction .3 \
                        --no-ewens-priors \
                        --use-mapping-quality \
                        --no-marginals \
                        --left-align-indels \
                        --use-mapping-quality \
                        --min-base-quality 25 \
                        --min-supporting-quality 30,30''',
                    region= region, **kwargs),
                vcf_outfile=output_prefix + 'vcf_compare.fb.eff.vcf',
                html_outfile=output_prefix + 'vcf_compare.fb.eff.html'))
    except StopIteration:
        return (region, [], None)

    records = []
    
    for record in reader:
        records.append(record)
    
    return (region, records, 
        {'metadata': reader.metadata, 'infos': reader.infos,
         'formats': reader.formats, 'filters': reader.filters,
         'samples': reader.samples})


def compare_vcf(bam_filenames, output_prefix, region, 
    genome= '/scratch/dbg/genomes/mg1655ref.fa', **kwargs):
    ''' compares many vcf callers in a region (to ease parallelization)
        by instantiating many vcf.reader objects and combining the output
        vcf files.
    '''
    
    
    out_stream = StringIO()
    SnpDef = namedtuple('SnpDef', ['CHROM', 'POS', 'REF', 'ALT'])
    
    snpdefs = defaultdict(dict)
    
    reader_list =get_reader_list(bam_filenames, output_prefix, region, 
        **kwargs)
    
    #for record, desc_set in descriptor_sets.items():
    #columns:
    # Position 'pos'
    # Effect 'eff'
    # Gene Name 'gene'
    # Depth: Ref/Alt 'depth'
    # seq_call: Ref/Het/Hom/-
    
    # ======================================================================
    # = Generate Columns for DataFrame ----------------------------------- =
    # ======================================================================
    columns = ['pos','ref','alt','nhom','nhet','genes','eff','gene','depth']
    
    # add reader-specific columns
    for reader_name, reader in reader_list.items():
    
        #add sample and reader name to columns
        columns.extend(['_'.join([reader_name,sample]) 
            for sample in reader.samples])
    
        #add records to snpdef dict
        for record in reader:
            snpdef = SnpDef(record.CHROM, record.POS, record.REF, 
                frozenset(record.ALT))
            snpdefs[snpdef][reader_name] = record
            #print record, ' added to '+reader_name
    
    print >> out_stream, '\t'.join(columns)
    
    #make a list to hold records
    records = []
    
    # ======================================================================
    # = FOR EACH RECORD -------------------------------------------------- =
    # ======================================================================
    for snpdef, readers in sorted(snpdefs.items(), key=lambda i: i[0].POS):
        first_rec = readers.values()[0]
        
        #FILTERING STEP:
        if not ('HIGH' in first_rec.INFO['EFF']
          or 'MODERATE' in first_rec.INFO['EFF']):
            #continue
            pass
        
        output = {}
        counts = defaultdict(int)
    
        output['pos'] = str(first_rec.POS)
        output['ref'] = str(first_rec.REF)
        output['alt'] = str(first_rec.ALT)
        output['eff'] = ''
        output['genes'] = ''
    
        effects = first_rec.INFO['EFF'].split(',')
        for eff in effects:
            eff = eff.replace('(','|').split('|')
            if eff[0] == 'INTRON': 
                #output['genes'] += eff[5]+','
                continue
            else:
                output['eff'] += (':'.join(eff[0:5])+';')
                output['genes'] += eff[5]+','
    
        if 'mp' in readers:
            record = readers['mp']
            output['depth'] = str(sum(record.INFO['DP4'][0:1])) + '/' + str(
                sum(record.INFO['DP4'][2:3]))
        else:
            record = readers[readers.keys()[0]]
            output['depth'] = str(record.INFO['RO']) + '/' + str(
                record.INFO['AO'])
    
        for reader_name in readers:
            for call in readers[reader_name].samples:
                sample = call.sample
                column = '_'.join([reader_name, sample])
            
                if call.is_variant:
                    if call.is_het: 
                        output[column] = 'Het'
                        counts['nhet'] += 1
                    else: 
                        output[column] = 'Hom'
                        counts['nhom'] += 1
                else: output[column] = 'Ref'
        
        
        for count,num in counts.items():
            output[count] = str(num)
        
        print >> out_stream, '\t'.join(
            output.setdefault(col, '-') for col in columns)
        
        #populate row in record dict
        records.append(output)
    
    if len(records) > 0:
        df = DataFrame(records, columns= columns)
        df = df.fillna('-')
        return (region, df, out_stream.getvalue())
    else:
        return (region, None, out_stream.getvalue())


def get_reader_list(bam_filenames, output_prefix, region, **kwargs):
    
    fn_string = ' '.join(bam_filenames)
    
    kwargs['region'] = region
    
    reader_list = {}
    
    try:
        reader_list['fb'] = vcf.Reader(
            call_snpeff(
                call_freebayes(fn_string,
                    fb_params=''' \
                        --pvar 0.001 \
                        --ploidy 2 \
                        --min-alternate-fraction .3 \
                        --no-ewens-priors \
                        --use-mapping-quality \
                        --no-marginals \
                        --left-align-indels \
                        --use-mapping-quality \
                        --min-base-quality 25 \
                        --min-supporting-quality 30,30''',
                    **kwargs),
                vcf_outfile=output_prefix + 'vcf_compare.fb.eff.vcf',
                html_outfile=output_prefix + 'vcf_compare.fb.eff.html'))
    except StopIteration:
        pass
    
    # try:
    #     reader_list['fbl'] = vcf.Reader(
    #         call_snpeff(
    #             call_freebayes(fn_string, 
    #                 fb_params=''' \
    #                     --pvar 0.001 \
    #                     --ploidy 2 \
    #                     --min-alternate-count 4 \
    #                     --no-ewens-priors \
    #                     --no-marginals''',
    #                 **kwargs),
    #             vcf_outfile=output_prefix + 'vcf_compare.fbl.eff.vcf',
    #             html_outfile=output_prefix + 'vcf_compare.fbl.eff.html'))
    # except StopIteration:
    #     pass
    # 
    # try:
    #     reader_list['mp'] = vcf.Reader(
    #         call_mpileup(call_mpileup(fn_string, **kwargs),
    #             vcf_outfile=output_prefix + '.vcf_compare.mp.eff.vcf',
    #             html_outfile=output_prefix + '.vcf_compare.mp.eff.html')) 
    # except StopIteration:
    #     pass
    
    return reader_list       



def make_vcf_df(output_prefix):
    
    #THIS IS BROKEN NOW THAT I HAVE UPDATED PyVCF =(
    #from prev run:
    #output_prefix = '/scratch/dbg/recoli/notes/'
    #df = pickle.load(open(output_prefix + 'data_frame.p', 'rb'))
   
    amberfile = open(output_prefix + 'ambers.sorted.sites')
    ambers = []
    for line in amberfile: ambers.append(int(line))
    ambers = Series(ambers)
    df['AMBER'] = df.POS.isin(ambers)
    df['MQM'] = [r.INFO['MQM'][0] for r in df.RECORD]
    df['PAIRED'] = [r.INFO['PAIRED'][0] for r in df.RECORD]
    df['DPRA'] = [r.INFO['PAIRED'][0] for r in df.RECORD]
    
    #EFF FIELD
    def get_eff(record):
        
        genes = set()
        eff_sev = set()
        eff_type = set()
        mut_pos = set()
        
        if 'EFF' in  record.INFO:
            effects = record.INFO['EFF'].split(',')
            for eff in effects:
                eff = eff.replace('(','|').split('|')
                if eff[0] == 'INTRON': 
                    continue
                eff_type.add(eff[0])
                eff_sev.add(eff[1])
                genes.add(eff[5])
                mut_pos.add(eff[4])
        
        return ['|'.join(i) for i in (genes, eff_sev, eff_type, mut_pos)]
    
    eff_fields = [get_eff(record) for record in df.RECORD]
    (genes, eff_sev, eff_type, mut_pos) = zip(*eff_fields)
    df['EFF_GENE'] = genes
    df['EFF_SEV'] = eff_sev
    df['EFF_TYPE'] = eff_type
    df['MUT_POS'] = mut_pos
    
    #mark parents and children
    df['STRAIN'] = [sample.split('_')[2] for sample in df.SAMPLE]
    
    parents = defaultdict(set)
    children = defaultdict(set)
    strain_num = {}
    

    samples = np.unique(df.SAMPLE)
    strains = np.unique(df.STRAIN)
    
    #get parent and child matching strains
    for strain in np.unique(df.STRAIN):
        
        match = re.match(r'(c|rec)(\d+)',strain)
        if not match:
            if strain == 'mg1655':
                strain_num[strain] = -1
                children[strain].add('ecnr2')
            if strain == 'ecnr2':
                strain_num[strain] = 0
                [children[strain].add('rec'+'%02d' % (i,)) for i in range(1,31)]
            continue
        
        (typ,num) = re.match(r'(c|rec)(\d+)',strain).groups()
        num = int(num)
        if typ == 'rec':
            children[strain].add('c' + '%02d' % ((int(num) + 1) / 2))
            strain_num[strain] = int(num)
        elif typ == 'c' and num < 31:
            children[strain].add('c' + '%02d' % ((int(num) + 33) / 2))
            strain_num[strain] = int(num) + 32
        elif typ == 'c' and num == 31:
            strain_num[strain] = 63
            pass #no children
        for child in children[strain]:
            parents[child].add(strain)
    
    #identify parent/child samples as well as strains
    strain_samples = dict([(strain, 
        list(samples[['_'+strain in sample for sample in samples]])) 
            for strain in strains])
    
    strain_parents = defaultdict(set)
    strain_children = defaultdict(set)
    
    for strain in strains:
        for parent in parents[strain]:
            for parent_sample in strain_samples[parent]:
                for sample in strain_samples[strain]:
                    strain_parents[sample].add(parent_sample)
            
        for child in children[strain]:
            for child_sample in strain_samples[child]:
                for sample in strain_samples[strain]:
                    strain_children[sample].add(child_sample)
    
    df_var = df[df.VAR != False]
    
    in_child_col = []
    in_parent_col = []
    
    for (sample, var, record) in zip(df_var.SAMPLE, df_var.VAR, df_var.RECORD):
        
        in_parents = 0
        for parent in strain_parents[sample]:
            if record.genotype(parent).is_variant != False:
                in_parents += 1
        
        in_children = 0
        for child in strain_children[sample]:
            if record.genotype(child).is_variant != False:
                in_children += 1
        
        in_parent_col.append(in_parents)
        in_child_col.append(in_children)
    
    df_var['IN_PARENT'] = in_parent_col
    df_var['IN_CHILD'] = in_child_col
    df_var['STRAIN_NUM'] = [strain_num[strain] for strain in df_var.STRAIN]
    
        
    df_var_sm = df_var[['CHROM', 'POS', 'SAMPLE', 'REF', 'TYPE', 'CALL', 
        'VAR', 'HET', 'BASE','AMBER', 'MQM', 'PAIRED', 'DPRA', 'EFF_GENE', 
        'EFF_SEV', 'EFF_TYPE', 'MUT_POS', 'STRAIN', 'STRAIN_NUM',
        'IN_PARENT', 'IN_CHILD']]
    df_var_sm.to_csv('all_snp_data.csv')
        
def make_pindel_VCF(output_prefix):
    #LOAD SVs from PINDEL=====================================================
    pindel_dict_list = []
    pindel_files = [output_prefix +  '../pindel/recoli_pindel_EVENTS.vcf',
        output_prefix +  '../pindel/recoli_pindel_BP.vcf']
    
    for pindel_file in pindel_files:
        
        out_stream = StringIO()
        pf_reader = vcf.Reader(call_snpeff(open(pindel_file), file=True))
        writer = vcf.Writer(out_stream, pf_reader)
        
        for record in pf_reader:
            writer.write_record(record)
            last_newline = out_stream.getvalue()[:-2].rfind('\n')
            record_line = out_stream.getvalue()[last_newline+1:-2] + '\n'

            gts = list(record.REF) + record.ALT

            for sample in pf_reader.samples:

                sample_call = record.genotype(sample)

                if sample_call.data.GT == None:

                    dict_record = {
                        'CHROM': record.CHROM, 
                        'POS': record.POS,
                        'SAMPLE': sample,
                        'REF': record.REF,
                        'TYPE': record.var_subtype,
                        'CALL': None,
                        'VCFLINE': record_line,
                        'RECORD': record,
                        'SVTYPE': record.INFO['SVTYPE'],
                        'SVLEN': record.INFO.get('SVLEN', None),
                        'END': record.INFO.get('END',''),
                        'BPDIR': record.INFO.get('BPDIR','')
                    }

                else:    

                    dict_record = {
                        'CHROM': record.CHROM, 
                        'POS': record.POS,
                        'SAMPLE': sample,
                        'REF': record.REF,
                        'TYPE': record.var_subtype,
                        #'CALL': gt_num,
                        'VAR': sample_call.is_variant,
                        #'HET': '', #no het call for SVs
                        #'BASE': gts[gt_num],
                        'VCFLINE': record_line,
                        'RECORD': record,
                        'SVTYPE': record.INFO['SVTYPE'],
                        'SVLEN': record.INFO.get('SVLEN', None),
                        'END': record.INFO.get('END',''),
                        'BPDIR': record.INFO.get('BPDIR','')
                        
                    }
                
                if dict_record['SVLEN'] in (None,0) or abs(dict_record['SVLEN']) > 3:
                    pindel_dict_list.append(dict_record)

    df = DataFrame(pindel_dict_list,
        columns=['CHROM','POS','END','SAMPLE','REF','TYPE','CALL','VAR',
                 'HET','BASE','VCFLINE','RECORD', 'SVLEN', 'BPDIR'])

    #EFF FIELD
    def get_eff(record):

        genes = set()
        eff_sev = set()
        eff_type = set()
        mut_pos = set()

        if 'EFF' in  record.INFO:
            effects = record.INFO['EFF'].split(',')
            for eff in effects:
                eff = eff.replace('(','|').split('|')
                if eff[0] == 'INTRON': 
                    continue
                eff_type.add(eff[0])
                eff_sev.add(eff[1])
                genes.add(eff[5])
                mut_pos.add(eff[4])

        return ['|'.join(i) for i in (genes, eff_sev, eff_type, mut_pos)]

    eff_fields = [get_eff(record) for record in df.RECORD]
    (genes, eff_sev, eff_type, mut_pos) = zip(*eff_fields)
    df['EFF_GENE'] = genes
    df['EFF_SEV'] = eff_sev
    df['EFF_TYPE'] = eff_type
    df['MUT_POS'] = mut_pos

    #mark parents and children
    df['STRAIN'] = [sample.split('_')[2] for sample in df.SAMPLE]

    parents = defaultdict(set)
    children = defaultdict(set)
    strain_num = {}


    samples = np.unique(df.SAMPLE)
    strains = np.unique(df.STRAIN)

    #get parent and child matching strains
    for strain in np.unique(df.STRAIN):

        match = re.match(r'(c|rec)(\d+)',strain)
        if not match:
            if strain == 'mg1655':
                strain_num[strain] = -1
                children[strain].add('ecnr2')
            if strain == 'ecnr2':
                strain_num[strain] = 0
                [children[strain].add('rec'+'%02d' % (i,)) for i in range(1,31)]
            continue

        (typ,num) = re.match(r'(c|rec)(\d+)',strain).groups()
        num = int(num)
        if typ == 'rec':
            children[strain].add('c' + '%02d' % ((int(num) + 1) / 2))
            strain_num[strain] = int(num)
        elif typ == 'c' and num < 31:
            children[strain].add('c' + '%02d' % ((int(num) + 33) / 2))
            strain_num[strain] = int(num) + 32
        elif typ == 'c' and num == 31:
            strain_num[strain] = 63
            pass #no children
        for child in children[strain]:
            parents[child].add(strain)

    #identify parent/child samples as well as strains
    strain_samples = dict([(strain, 
        list(samples[['_'+strain in sample for sample in samples]])) 
            for strain in strains])

    strain_parents = defaultdict(set)
    strain_children = defaultdict(set)

    for strain in strains:
        for parent in parents[strain]:
            for parent_sample in strain_samples[parent]:
                for sample in strain_samples[strain]:
                    strain_parents[sample].add(parent_sample)

        for child in children[strain]:
            for child_sample in strain_samples[child]:
                for sample in strain_samples[strain]:
                    strain_children[sample].add(child_sample)

    df_var = df[(df.VAR == True) & (df.TYPE != 'RPL')]

    in_child_col = []
    in_parent_col = []

    for (sample, var, record) in zip(df_var.SAMPLE, df_var.VAR, df_var.RECORD):

        in_parents = 0
        for parent in strain_parents[sample]:
            if record.genotype(parent).is_variant != False:
                in_parents += 1

        in_children = 0
        for child in strain_children[sample]:
            if record.genotype(child).is_variant != False:
                in_children += 1

        in_parent_col.append(in_parents)
        in_child_col.append(in_children)

    df_var['IN_PARENT'] = in_parent_col
    df_var['IN_CHILD'] = in_child_col
    df_var['STRAIN_NUM'] = [strain_num[strain] for strain in df_var.STRAIN]

    df_var_sm = df_var[['CHROM', 'POS', 'END', 'SAMPLE', 'REF', 'TYPE', 
        'VAR', 'HET', 'BASE', 'EFF_GENE', 'BPDIR',
        'EFF_SEV', 'EFF_TYPE', 'MUT_POS', 'STRAIN', 'STRAIN_NUM',
        'IN_PARENT', 'IN_CHILD']]
    df_var_sm.to_csv('pindel_sv_data.csv')


def compare_vcf_parallel(bam_filenames, output_prefix, 
    genome= '/scratch/dbg/genomes/mg1655ref.fa', **kwargs):
    
    #genome_size = 4639675
    genome_size = 60000
    chunk_size = 5000
    
    pool = Pool(12)
    jobs = {} #jobs
    
    chunk_outputs = []
    
    chunks = zip(range(0, (genome_size-chunk_size), chunk_size), 
        range(chunk_size+1,genome_size,chunk_size))
    
    for region in chunks:
        
        #kwargs['region'] = region
                    
        jobs[region] = pool.apply_async(compare_vcf_calls,
            (sys.argv[1:], output_prefix, region,), 
            kwargs, callback= chunk_outputs.append)
        # chunk_outputs.append(
        #     compare_vcf_calls(sys.argv[1:], output_prefix, region, **kwargs))
    
    while True:
    
        j_done = sum([j.ready() for j in jobs.values()])
        j_tot = len(chunks)
        
        print >> sys.stderr, '\t%2d / %2d jobs completed.\r' % (
                j_done, j_tot),
        sys.stderr.flush()
    
        if j_done == j_tot:
            for job in jobs.values():
                try: assert job.successful() 
                except AssertionError: job.get()
            print >> sys.stderr, '\n'
            break
            
        else:
            time.sleep(1)
    
    (regions, dfs, plaintxts) = zip(*chunk_outputs)
    df = concat(dfs, ignore_index=True)
    
    ## find records where lax freebayes and mp find something, but conserv
    ## freebayes does not
    bad_positions = df[( (df['fb_recoli_lib3_c31_322A'] == '-') & ((df['fbl_recoli_lib3_c31_322A'] != '-') | (df['mp_recoli_lib3_c31_322A'] != '-') ))]['pos']
    
    ## remove positions where strong freebayes has a snp called but for a 
    ## different record
    bad_positions = np.unique(bad_positions[np.logical_not(
        bad_positions.isin(
            df[(df['fb_recoli_lib3_c31_322A'] != '-')]['pos']))]) 
    
    for pos in bad_positions:
        print 'Chromosome\t%d' % pos    


if __name__ == "__main__":

    output_prefix = '/scratch/dbg/recoli/notes/'
    kwargs = {}
    #kwargs['bam_filenames'] = sys.argv[1:]
    kwargs['fasta_reference_filename'] = '/scratch/dbg/genomes/mg1655ref.fa'
    
    #write_vcf_parallel(sys.argv[1:], output_prefix, **kwargs)
    #make_vcf_df(output_prefix)
    make_pindel_VCF(output_prefix)

# vcfttview \
#     -B $proj_dir/align/recoli_lib3_c31_322B/recoli_lib3_c31_322B.realign.bam \
#     -R /scratch/dbg/genomes/mg1655ref.fa \
#     <( perl -pe "s/\s+[0]+/\t/; s/ +/\t/g;" /tmp/badpos.vcf)

# reader = vcf.Reader(open('/scratch/dbg/recoli/notes/vcf_compare.fb.eff.vcf'))
# writer = vcf.Writer(out_stream, reader)
    