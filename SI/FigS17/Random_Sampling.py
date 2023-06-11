import numpy as np
import multiprocessing
import pandas as pd
from scipy.integrate import odeint
from Latin_Hypercube import LHSample
from Equation import func
from Type_Characterization2 import Type_Characterization,WithinPhys
import warnings
import time

warnings.filterwarnings('error')

def get_mode(Para):
    Time=np.arange(0, 80, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    try:
        result=odeint(func, initial, Time, args=(Para,))
    except:
        return 100, Para
    wp=WithinPhys(result)
    if wp:
        mode=Type_Characterization(result)
        return mode,Para
    else:
        #print(constraints)
        return 100,Para

if __name__=='__main__':
    sample_indices1 = np.arange(0, 25, 1)
    sample_indices1 = np.delete(sample_indices1, [8, 9, 11, 12, 13, 14, 15, 17, 18, 19, 21, 22, 23, 10, 20, 24])
    sample_indices2 = np.hstack((np.array([8, 9, 11, 12, 13, 14, 15, 17, 18, 19, 21, 22, 23]), [30, 34, 35],
                                 np.arange(36, 52, 1)))
    sample_indices3 = [150, 159, 160] + list(np.arange(128, 132, 1)) + list(np.arange(27, 30, 1)) + list(
        np.arange(31, 34, 1))
    para = np.loadtxt('Para.txt')
    n = len(para)
    LogParaBound = []
    for i in range(len(para)):
        a = np.log10(para[i])
        if i in sample_indices1:
            LogParaBound.append([a - np.log10(5), a + np.log10(5)])
        elif i in sample_indices2:
            LogParaBound.append([a - np.log10(2), a + np.log10(2)])
        elif i in sample_indices3:
            LogParaBound.append([a - np.log10(5), a + np.log10(5)])
        elif i not in [25,61,62,151,52]:
            LogParaBound.append([a - np.log10(2), a + np.log10(2)])
        else:
            LogParaBound.append([a, a])
    LogParaBound = np.array(LogParaBound)

    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-1))
    N=2*10**3
    time0=time.perf_counter()
    LogPara_matrix=np.array(LHSample(len(LogParaBound), LogParaBound, N))
    Para_matrix=np.power(10,LogPara_matrix)
    Para_mode = [[] for i in range(5)]
    count=0
    df=[]
    for y in pool.imap(get_mode, Para_matrix):
        # print
        ratio=count/N
        rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
        rat_str=''.join(rat_str)
        print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
        count+=1
        mode=y[0]
        Para=y[1]
        if mode<=4:
            Para_mode[mode].append(Para)
        if 1<=mode<=4:
            Para_mode[0].append(Para)
        if count%1000==0:
            print('\n processes completed:',count,'Within_Phys_Range:',len(Para_mode[0]),
                  'Mode 1:%i/%i'%(len(Para_mode[1]),len(Para_mode[0])),
                  'Mode 2:%i/%i' % (len(Para_mode[2]), len(Para_mode[0])),
                  'Mode 3:%i/%i'%(len(Para_mode[3]),len(Para_mode[0])),
                  'Mode 4:%i/%i' % (len(Para_mode[4]), len(Para_mode[0])),
                  'time used:%1.1f s'%(time.perf_counter()-time0))
    for mode in range(0, 5, 1):
        np.savetxt('rSpl_Para%i.txt'%mode, Para_mode[mode])
    pool.close()