import subprocess
import tempfile
import util

#count_file1 = iter('1 a\n1 b\n3 c\n4 d\n1 q'.split('\n'))
#count_file2 = iter('1 a\n2 b\n3 c\n4 q\n5 z'.split('\n'))
 

def ngs_mcount_reduce(count_file1, count_file2, **kwargs):
    
    pass
    
    #NOTE: THIS IS CURRENTLY BROKEN!!
    count_file1 = iter(count_file1)
    count_file2 = iter(count_file2)
    
    total_merged_count = 0
    total_count_1 = 0
    total_count_2 = 0
    
    line1 = count_file1.next()
    line2 = count_file2.next()
    
    (count1, seq1) = line1.split()
    (count2, seq2) = line2.split()

    while True:
        try:
            cmpval = cmp(seq1, seq2)
    
            if cmpval == 0:
                merged_count = int(count1) + int(count2)
                print ' '.join((str(merged_count), seq1))
                
                line1 = count_file1.next()
                (count1, seq1) = line1.split()
                
                line2 = count_file2.next()
                (count2, seq2) = line2.split()
                
                total_count_1 += count_1
                total_count_2 += count_2
                total_merged_count += (int(count_1) + int(count_2))
                continue
    
            elif cmpval < 0 :
                print line1.lstrip().rstrip()
                line1 = count_file1.next()
                (count1, seq1) = line1.split()
                
                total_count_1 += count_1
                total_merged_count += int(count_1)
                continue
        
            elif cmpval > 0:
                print line2.lstrip().rstrip()
                line2 = count_file2.next()
                (count2, seq2) = line2.split()
                
                total_count_2 += count_2
                total_merged_count += int(count_2)
                continue
            
        except StopIteration:
            #one of the files has ended before the other,
            # write the remainder of both files to the merged one
            
            remaining1 = map(lambda line: (line, line.split()), 
                list(count_file1))
            remaining2 = map(lambda line: (line, line.split()), 
                list(count_file2))
                
            #should have three values, line, count, seq
            for vals in remaining1, remaining2:
                assert len(vals) == 3
            
            if cmpval > 0:
                print line1.lstrip().rstrip()
                total_count_1 += count_1
                total_merged_count += int(count_1)
                
                for line, count, seq in remaining1:
                    print line1.lstrip().rstrip()
                    total_merged_count += count
                    total_count_1 += count
                
                for line, count, seq in remaining2:
                    total_merged_count += count
                    total_count_1 += count
                    
                
            elif cmpval < 0:
                print '\n'.join([line2.lstrip().rstrip()] + remaining2)
                total_count_1 += count_2
                total_merged_count += int(count_2)
                total_count_2 += count_2
                total_merged_count += int(count_2)
                
                
            elif cmpval == 0:
                if remaining1: print '\n'.join(remaining1)
                if remaining2: print '\n'.join(remaining2)
            break

def ngs_mcount_map(seq_file, **kwargs):
    ''' for multi-thread counting, use the shell to count unique reads per
        tile/bin, then return a file-like obj with the counts for that bin.'''

    count_file = tempfile.SpooledTemporaryFile(max_size=8e5)

    command = 'grep -Pi \'^[ATGC]+$\' | sort | uniq -c'
    
    child = subprocess.call(command,
        stdout=count_file,
        stdin=seq_file,
        shell=True,
        executable='/bin/bash')
    
    count_file.seek(0)
    
    return count_file
    
if __name__ == '__main__':
    #fn1 = '/Volumes/gmc/scratch/dbg/ecre/ct/202_bins/202_bins.1.001.MT.seq.gz'
    #fn2 = '/Volumes/gmc/scratch/dbg/ecre/ct/202_bins/202_bins.1.002.MT.seq.gz'
    fn1 = '/scratch/dbg/ecre/ct/202_bins/202_bins.1.001.MT.seq.gz'
    fn2 = '/scratch/dbg/ecre/ct/202_bins/202_bins.1.002.MT.seq.gz'
    
    file1 = util.open_zcat(fn1)
    file2 = util.open_zcat(fn2)
        
    reduce(ngs_mcount_reduce,map(ngs_mcount_map,(file1, file2)))