import numpy as np
import multiprocessing
import pandas as pd
from scipy.integrate import odeint
from Latin_Hypercube import LHSample
from Equation import func
from Type_Characterization2 import Type_Characterization_test,WithinPhys
import warnings
import time

warnings.filterwarnings('error')

def get_mode(Para,a,b):
    Time=np.arange(0, 80, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    result=odeint(func, initial, Time, args=(Para,))
    wp,constraints=WithinPhys(result)
    if wp:
        mode=Type_Characterization_test(result,a,b)
        return mode,Para,constraints
    else:
        #print(constraints)
        return 100,Para,constraints

def get_mode_merge(args):
    return get_mode(*args)

if __name__=='__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-1))
    ab=[[800,1600],[1200,2400]]
    for ii in range(2):
        a,b=ab[ii]
        N=2*10**3
        time0=time.perf_counter()
        LogParaBound=np.loadtxt('logparabound.txt')
        LogPara_matrix=np.array(LHSample(len(LogParaBound), LogParaBound, N))
        Para_matrix=np.power(10,LogPara_matrix)
        args=[[Para_matrix[i,:],a,b] for i in range(N)]
        Para_mode = [[] for i in range(5)]
        count=0
        for y in pool.imap(get_mode_merge, args):
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
            if count%500==0:
                print('\n processes completed:',count,'Within_Phys_Range:',len(Para_mode[0]),
                      'Mode 1:%i/%i'%(len(Para_mode[1]),len(Para_mode[0])),
                      'Mode 2:%i/%i' % (len(Para_mode[2]), len(Para_mode[0])),
                      'Mode 3:%i/%i'%(len(Para_mode[3]),len(Para_mode[0])),
                      'Mode 4:%i/%i' % (len(Para_mode[4]), len(Para_mode[0])),
                      'time used:%1.1f s'%(time.perf_counter()-time0))
        for mode in range(0, 5, 1):
            np.savetxt('./Test_Criteria/Spl_Para%i_%i.txt'%(ii,mode), Para_mode[mode])
    pool.close()