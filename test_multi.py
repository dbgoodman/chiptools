from multiprocessing import Process, Manager, Pool, Lock
import pprint

def f(seqs_in_bin, bin, bins, seq_dict, bin_list, seq_dict_lock, bin_locks,
      bin_list_lock):
    
    bin_idx = bins.index(bin)
    
    for seq in seqs_in_bin:
        bin_locks[bin].acquire()
        #seq_dict_lock.acquire()
        #bin_list_lock.acquire()
        
        
        if seq in seq_dict:
            bin_list_idx = seq_dict[seq][0] + bin_idx
            #bin_locks[bin].acquire()
            bin_list[bin_list_idx] += 1
            #bin_locks[bin].release()
        else:
            #seq_dict_lock.acquire()
            bin_list_lock.acquire()
            new_bin = len(bin_list)
            new_bins = [0] * len(bins)
            bin_list.extend(new_bins)
            #bin_list_lock.release()
            seq_dict[seq] = (new_bin, new_bin+len(bins))
            #seq_dict_lock.release()
            #bin_locks[bin].acquire()
            bin_list[new_bin+bin_idx] += 1
            #bin_locks[bin].release()
        
        bin_locks[bin].release()
        #seq_dict_lock.release()
        #bin_list_lock.release()
        
            

if __name__ == '__main__':
    mgr = Manager()
    
    bins = ['a','b','c','d','e']
    num_bins = len(bins) 
    seq_dict = mgr.dict() #seq, bin tuple to idx number
    bin_list = mgr.list() #idx to counts
    
    bin_locks = mgr.dict() #lock for each bin
    seq_dict_lock = mgr.RLock() #lock for seq_dict
    bin_list_lock = mgr.RLock() #lock for whole bin list
    
    for bin in bins:
        bin_locks[bin] = mgr.RLock()
    
    seqs = mgr.dict()
    
    seqs['a'] = ['GATC','GATC','ATC','AAG']
    seqs['b'] = ['GATC','CCTC','TTC','GGT']
    seqs['c'] = ['GATC','CCTC','GCG','TCT']
    seqs['d'] = ['GATC','AAGT','CTC','CTC']
    seqs['e'] = ['GATC','AAGT','CCT','CCT']
    
    ppool = list()
    
    for bin, seqs_in_bin in seqs.items():
        print bin
        proc = Process(target=f, args=( 
            seqs_in_bin, bin, bins, seq_dict, bin_list, 
            seq_dict_lock, bin_locks, bin_list_lock))
        
        proc.start()
        ppool.append(proc)
    
    for proc in ppool:
        proc.join(timeout=1)
    
    print seq_dict
    print bin_list

# from multiprocessing import Process, Manager
# def f(m, d):
#     if not 'f' in d:
#         d['f'] = m.dict()
#     else:
#         d['f'] = 1
# 
# m = Manager()
# d = m.dict()
# 
# pool = list()
# 
# for i in range(5):
#     p = Process(target=f, args=(m,d))
#     pool.append(p)
#     p.start()
#     
# for p in pool:
#     p.join()
#     
# print d
# print d['f']