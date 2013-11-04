#http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic

import multiprocessing
import multiprocessing.pool
import math


def parallel_reduce(func, iterable, processes= 4, args=(), kwargs={}):
    #iterable = range(1,11)
    comp_stack = list(iterable)
    pair_list = []
    
    pool = multiprocessing.pool.Pool(processes)
    #print comp_stack

    while len(comp_stack) > 1:
        while len(comp_stack) > 1:
            pair_list.append((comp_stack.pop(), comp_stack.pop()))
            
        #print pair_list
    
        results = []
        while len(pair_list) > 0:    
            pair = pair_list.pop()
            results.append(pool.apply_async(func, [pair]))
    
        #print results 
                
        while True:
            if all([result.ready() for result in results]): break
        
        comp_stack = [result.get() for result in results]
    
    return comp_stack

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class RecursivePool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


def test_recursive_parallel_reduce(workers = 5):
    
    pool = RecursivePool()
    
    ranges = [range(1, 5), range(2, 9), range(3, 7)]
    
    print ranges
    
    results = []

    for myrange in ranges:
        pool.apply_async(parallel_reduce, [sum, myrange], 
            callback= results.append)

    pool.close()
    pool.join()

    print results

if __name__ == '__main__':
    test_recursive_parallel_reduce()