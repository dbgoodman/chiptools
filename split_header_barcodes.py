from collections import defaultdict
import os
import gzip
import sys
import time

samples = {'ACCTGAAT':'%s/202_hsrna/s_G1_L001_R%d_%02d.fastq.1.gz',
           'CAAGGCAT':'%s/202_hsrna/s_G1_L001_R%d_%02d.fastq.2.gz', 
           'GTTCATAT':'%s/203_hsrna/s_G1_L001_R%d_%02d.fastq.1.gz',
           'TATAAGAT':'%s/203_hsrna/s_G1_L001_R%d_%02d.fastq.2.gz'}

def getfilename(fn, fdir, read, split): 
    return fn % (fdir, read, split)

def split_reads(r1_in, r2_in, fdir, num_split, part):
    
    num_split = int(num_split)
    part = int(part) - 1
    
    #make input files
    r1_inf = gzip.open(r1_in)
    r2_inf = gzip.open(r2_in)
    
    missing_total = 0
    read_totals = defaultdict(int)
    
    #make output files
    fobjs = {}
    for bc, fn in samples.items():
        for split in range(num_split):
            for read in (1,2):
                full_fn = getfilename(fn, fdir, read, 
                    (part*num_split) + split + 1)
                print "Making file "+full_fn
                fobjs[(bc,split,read)] = gzip.open(full_fn, 'w')
    
    split = 0
    try:
        while 1:
            r1 = [r1_inf.next(),r1_inf.next(),r1_inf.next(),r1_inf.next()]
            r2 = [r2_inf.next(),r2_inf.next(),r2_inf.next(),r2_inf.next()]
            
            if r1[0][:-11] != r2[0][:-11]:
                raise ValueError("HEADER MISMATCH: %s vs. %s" % (r1[0], r1[0])) 
            barcode = r1[0][-11:-3]
            
            if barcode not in samples.keys():
                missing_total += 1
            else:
                read_totals[barcode] += 1
                fobjs[(barcode,split,1)].write(''.join(r1))
                fobjs[(barcode,split,2)].write(''.join(r2))
            
                #cycle through split files
                split += 1
                if split >= num_split: split = 0
            
    except StopIteration:
        pass
    
    total_read_count = sum(read_totals.values()) + missing_total
    print "Sample\tBarcode\tCount\tPercent"
    for barcode in sorted(samples.keys()):
        print "%s\t%s\t%s\t%.2f%%" % (samples[barcode], barcode, 
            read_totals[barcode], 
            (float(read_totals[barcode])/float(total_read_count))*100)
    print "%s\t%s\t%s\t%.2f%%" % ('unmatched', None, missing_total, 
        (float(missing_total)/float(total_read_count))*100)

if __name__ == '__main__':
    sys.exit(split_reads(sys.argv[1], sys.argv[2], sys.argv[3], 
        sys.argv[4], sys.argv[5]))