import argparse
import sys

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

import fasta
import align
import simchip
import ngs2

RESTR1 = 'GGCGCGCC'
RESTR2 = 'CATATG'

def main():
    # parse user input
    args = parser.parse_args()
    args.func(**args.__dict__)

# ============================================================================
# = TOP LEVEL PARSER FOR CHIPTOOLS
# ============================================================================
parser = argparse.ArgumentParser(prog='chiptools')
parser.add_argument('--verbose', '-v', action='count')

# ============================================================================
# = SUBPARSERS FOR CHIPTOOLS
# ============================================================================
subparsers = parser.add_subparsers(help='available subcommands')
#EXTRACT_SANGER---------------------------------------------------------------
parser_extract_sanger = subparsers.add_parser(
    'extract_sanger',
    help='''
    Extracts Sanger Sequences from a Directory of *.seq files into a
    single FASTA file. Keeps only the subset that contains a sequence
    between two specified restriction sites''')
        
parser_extract_sanger.set_defaults(func=fasta.extract_sanger)
    
parser_extract_sanger.add_argument(
    '--regex', '-r', metavar='re', type=str, dest='regex_str',
    help='''regex for parsing plate (1)  and number (2), default to 
            r'r'^[\w-]*(\d+-\d*)-(\d+)' ''',
    default=r'^[\w-]*(\d+-\d*)-(\d+)')
    
parser_extract_sanger.add_argument(
    '--output', '-o', type=argparse.FileType('w'), 
    help='''FASTA file output''',
    default=sys.stdout)

parser_extract_sanger.add_argument(
    '--rs1', metavar='restr_site_1', type=str, 
    help='''first restriction site, defaults to '''+RESTR1,
    default=RESTR1)

parser_extract_sanger.add_argument(
    '--rs2', metavar='restr_site_2', type=str, 
    help='''second restriction site, defaults to '''+RESTR2,
    default=RESTR2)

parser_extract_sanger.add_argument(
    'seq_path', metavar='seq_file_path',
    help='''directory containing seq files''')

parser_extract_sanger.add_argument(
    '--revcom','-R', action='store_true',
    dest='revcom', 
    help='''if flagged, reverse-complement the sanger reads after trimming''')

#EXTRACT_LIBRARY--------------------------------------------------------------
parser_extract_library = subparsers.add_parser(
    'extract_library',
    help='''
    Extract Library Members from a FASTA File containing reference
    sequences, trimming them between two specified restriction sites''')

parser_extract_library.set_defaults(func=fasta.extract_library)

parser_extract_library.add_argument(
    '--output', '-o', metavar='output_FASTA', 
    type=argparse.FileType('w'), 
    help='''output FASTA file''',
    default=sys.stdout)

parser_extract_library.add_argument(
    'input_filename', metavar='input_FASTA', 
    type=argparse.FileType('r'), 
    help='''library FASTA file''',
    default=sys.stdin)

parser_extract_library.add_argument(
    '--rs1', metavar='restr_site_1', type=str, 
    help='''first restriction site, defaults to '''+RESTR1,
    default=RESTR1)

parser_extract_library.add_argument(
    '--rs2', metavar='restr_site_2', type=str, 
    help='''second restriction site, defaults to '''+RESTR2,
    default=RESTR2)

parser_extract_library.add_argument(
    '--revcom','-R', action='store_true',
    dest='revcom', 
    help='''if flagged, reverse-complement the library after trimming''')

parser_extract_library.add_argument(
    '--regex', '-r', metavar='re', type=str, dest='re_str',
    help='''regex for parsing names, default to r'(.*+)' ''',
    default=r'(.*+)')

#PARSE_CLUSTER_STATS----------------------------------------------------------
parser_cluster_stats = subparsers.add_parser(
    'cluster_stats',
    help='''
    Print out a cluster summary files based on a cd-hit-2d cluster file.
    You can optionally include query and library fasta files (pre-trimmed)
    in order to get alignment, per-query and per-lib_member stats.''')

parser_cluster_stats.set_defaults(func=align.cluster_stats)

parser_cluster_stats.add_argument(
    '--output_prefix', metavar='output_prefix', 
    type=str, 
    help='''output file prefix. Defaults to ./output.*''',
    default='./output')

parser_cluster_stats.add_argument(
    '--query_seq', metavar='query_fasta_file', 
    type=argparse.FileType('r'), 
    help='''pre-trimmed query file''',
    default=None)

parser_cluster_stats.add_argument(
    '--lib_seq', metavar='library_fasta_file', 
    type=argparse.FileType('r'), 
    help='''pre-trimmed library file''',
    default=None)

parser_cluster_stats.add_argument(
    'cluster_file', metavar='cluster_file', 
    type=argparse.FileType('r'), 
    help='''FASTA file containing raw library to search against''',
    default=sys.stdin)

parser_cluster_stats.add_argument(
    '--local',
    action='store_true',
    help='''
    if flagged, use local alignment instead of global when doing 
    pairwise aligns''')    

parser_cluster_stats.add_argument(
     '--is_coding','-c', metavar='reading_frame', 
     dest='is_coding', type=int,
     help='''translate aligned NTs into codons in the aligns file. The reading
          frame must be positive.''')

#PARSE_cdhit2d----------------------------------------------------------
parser_cdhit2d = subparsers.add_parser(
    'cdhit2d',
    help='''
    Run cd-hit-2d on a lib file and a sanger seq dir, parse clusters, 
    perform pairwise alignment, and output stats to several files''')

parser_cdhit2d.set_defaults(func=align.cdhit2d)

parser_cdhit2d.add_argument(
    '--output', '-o', metavar='output_prefix', 
    type=str, dest='output_prefix',
    help='''output file prefix''',
    default='cdhit2d_output')

#lib_seqrecords, query_seqrecords, output_file
parser_cdhit2d.add_argument(
    'lib_seqrecords', metavar='library_file', 
    type=argparse.FileType('r'), 
    help='''FASTA file containing raw library to search against''')

parser_cdhit2d.add_argument(
    'query_seqrecords', metavar='seq_file_path',
    type=str,
    help='''FASTA file containing raw sanger sequences''')
    
parser_cdhit2d.add_argument(
    '--regex', '-r', metavar='re', type=str, 
    help='''regex for parsing plate (1)  and number (2), default to 
            r'r'^[\w-]*(\d+-\d*)-(\d+)' ''',
    default=r'^[\w-]*(\d+-\d*)-(\d+)')

parser_cdhit2d.add_argument(
    '--rs1', metavar='restr_site_1', type=str, 
    help='''first restriction site, defaults to '''+RESTR1,
    default=RESTR1)

parser_cdhit2d.add_argument(
    '--rs2', metavar='restr_site_2', type=str, 
    help='''second restriction site, defaults to '''+RESTR2,
    default=RESTR2)   

#SIMULATE CHIP----------------------------------------------------------------
parser_simchip = subparsers.add_parser(
    'simchip',
    help='''
    Given a list of fasta sequences, generate a fastq file that simulates the
    sequencing of a chip, including errors. The read depth will be normally 
    distributed among all of the sequences. Uses simNGS binary.''')

parser_simchip.set_defaults(func=simchip.simulate_chip)

parser_simchip.add_argument(
    '--num_fragments', '-f', metavar='#_fragments', 
    dest='num_fragments',
    type=int,
    help='''number of fragments to generate''',
    default=20000)

parser_simchip.add_argument(
    '--num_clones', '-c', metavar='#_clones', 
    dest='num_clones',
    type=int,
    help='''number of clones to generate (some will have mutations)''',
    default=2000)

parser_simchip.add_argument(
    'chip_fasta_fn', metavar='chip_fasta_file', 
    type=argparse.FileType('r'), 
    help='''FASTA file containing chip sequences''')

parser_simchip.add_argument(
    '--output', '-o', metavar='fasta_output_file', 
    type=argparse.FileType('w'), 
    help='''FASTA file with mutant clones from library 
            (defaults to stdout)''',
    default=sys.stdout)

#CLUSTER_NGS------------------------------------------------------------------
parser_cluster_ngs = subparsers.add_parser(
    'cluster_ngs',
    help='''
    Given raw paired end NGS reads, combine overlapping reads into merged
    reads and trim adapters and specified constant sequences. Then count
    the number of occurrences per unique read. Align the reads to the
    designed library, which must be pre-trimmed to match the NGS
    adapter sequences specified. 
    
    If multiple read1 and read2 files are given with -q, the counts for each
    are computed separately as well as together (i.e. for flow-seq bins).
    ''')
#FASTQ FILE MODES===============

parser_cluster_ngs.set_defaults(func=ngs2.ngs_cluster)

ngs_single_mode = parser_cluster_ngs.add_argument_group(
    title='1A: Single FASTQ file mode... ')

ngs_single_mode.add_argument(
    '--read1', '-1', metavar='read_1.fq.gz', 
    dest='fq_read1',
    type=argparse.FileType('r'),
    default= None,
    help='''The path to the gzipped FASTQ file containing R1 reads''')

ngs_single_mode.add_argument(
    '--read2', '-2', metavar='read_2.fq.gz', 
    dest='fq_read2',
    default= None,
    type=argparse.FileType('r'),
    help='''The path to the gzipped FASTQ file containing R2 reads''')

ngs_multi_mode = parser_cluster_ngs.add_argument_group(
    title='1B: ... or Multi FASTQ file mode')
    
ngs_multi_mode.add_argument(
    '--fastq', '-i', metavar='file.fq.gz', 
    dest='fq_files',
    type=argparse.FileType('r'),
    nargs='+',
    default= None,
    help='''The path to the gzipped FASTQ files containing all reads; parse
            using --regex, see help''')


#LIB FILE===============

ngs_library = parser_cluster_ngs.add_argument_group(
    title='2: Specify a library FASTA file')
    
ngs_library.add_argument(
    '--library', '-l', metavar='library_file',
    dest='lib_records', 
    required=True,    
    type=argparse.FileType('r'), 
    help='''The path to the FASTA file containing raw library to search 
        against (pre-trimmed)''')

#OPTIONS===============

ngs_options = parser_cluster_ngs.add_argument_group(
    title='3: Here are some optional settings')

ngs_options.add_argument(
    '--min_count', '-M', metavar='num_of_reads', 
    dest='tot_count',
    required=False,
    default=2,
    type=int,
    help='''only print out a unique read if it is seen at least this 
            many times. Default: 2''')

ngs_options.add_argument(
    '--min_length', '-L', metavar='min_length', 
    dest='minlen',
    required=False,
    type=int,
    help='''only keep a read if is at least N bases long after trimming. 
            Default: no minimum''')

ngs_options.add_argument(
    '--fwd_adapter', '-A', metavar='read1_adapter', 
    dest='adapter1',
    required=True,    
    type=str,
    help='''forward read primer/adapter sequence to trim as it would appear 
            at the end of a read''')

ngs_options.add_argument(
    '--rev_adapter', '-B', metavar='read2_adapter', 
    dest='adapter2',
    required=True,    
    type=str,
    help='''reverse read primer/adapter sequence to trim as it would appear 
            at the end of a read''')         

ngs_options.add_argument(
    '--revcom','-R', action='store_true',
    dest='revcom',
    help='''if flagged, reverse-complement the merged reads after SeqPrep,
            but before cd-hit-2d''')

ngs_options.add_argument(
    '--skip_finished','-s', action='store_true',
    dest='skip_finished',
    help='''if flagged, skip external steps that have completed in order to
            save time (useful for debugging only)''')

ngs_options.add_argument(
    '--regex', '-r', metavar='re', type=str, dest='regex_str',
    help='''If using multi FASTQ file mode, regex for parsing read and bin, 
            default:
            r'.*R(?P<read>[12])_(?P<tile>\w+).fastq.(?P<bin>\d+).gz' ''',
    default= r'.*R(?P<read>[12])_(?P<tile>\w+).fastq.(?P<bin>\d+).gz')

ngs_options.add_argument(
    '--rs1', metavar='restr_site_1', type=str, 
    help='''first restriction site, default: '''+RESTR1,
    default=RESTR1)

ngs_options.add_argument(
    '--rs2', metavar='restr_site_2', type=str, 
    help='''second restriction site, default: '''+RESTR2,
    default=RESTR2)

ngs_options.add_argument(
    '--manual_trim','-t', action='store_true',
    dest='no_manual_trim', 
    help='''do not trim at restriction sites after merging reads''')

ngs_options.add_argument(
    '--RNA', action='store_true',
    dest='RNA', 
    default=False,
    help='''UCLUST alignment and output for RNA instead of DNA''')

ngs_options.add_argument(
    '--no_uclust', action='store_true',
    dest='no_uclust',
    default=False,
    help='''Skip UCLUST step and only output unique fasta with bin counts''')

# parser_cluster_ngs.add_argument(
#     '--constant_seq','-c', metavar='constant_seq_to_trim', 
#     type=str, 
#     help='''a sequence that should be trimmed off of the read (from either
#             end) after merging and removing the adapters. Can be specified
#             multiple times.''')

#OUTPUT DIR===============
ngs_output_dir = parser_cluster_ngs.add_argument_group(
    title='4: Choose an output dir/prefix')
            
ngs_output_dir.add_argument(
    '--output', '-o', metavar='output_prefix', 
    dest='output_prefix', 
    type=str, 
    help='''
    path and file prefix for outputting each unique library member to various
    files.''',
    default=sys.stdout)
     
if __name__ == "__main__":
        main()
        
