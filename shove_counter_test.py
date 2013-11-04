import shove_counter
import shove
from multiprocessing import Process, Value, Array

def f(name, sc, sh, n, a):
    print 'hello', name
    sh['a'] = 5
    sc['b'] = 5
    print sc
    print sh
    n.value = 3.1415927
    for i in range(len(a)):
        a[i] = -a[i]

if __name__ == '__main__':
    sc = shove_counter.Counter()
    sh = shove.Shove(store='memory://', cache='memory://')
    num = Value('d', 0.0)
    arr = Array('i', range(10))
    
    p = Process(target=f, args=('sam', sc, sh, num, arr))
    p.start()
    p.join()
    print sc
    print sh
    print num.value
    print arr[:]


    