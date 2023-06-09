#Sampling the parameters in the range given by logparabound.py, run GetBound.py before running this one
#used with multiprocessing for parallel computation
#output txt files as 2D-array (n*len(Para)) for further analysis
#Para0: all the parameter sets within physiological range; Para1~3: mode 1~3; Para4: asymptomatic patients
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
    result=odeint(func, initial, Time, args=(Para,))
    wp,constraints=WithinPhys(result)
    if wp:
        mode=Type_Characterization(result)
        return mode,Para,constraints
    else:
        #print(constraints)
        return 100,Para,constraints

if __name__=='__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-1))
    N=2*10**3
    time0=time.perf_counter()
    LogParaBound=np.loadtxt('logparabound.txt')
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
        constraints=y[2]
        df.append(pd.DataFrame(constraints,index=[count]))
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
    df=pd.concat(df)
    df.loc['sum']=df.apply(lambda x: x.sum())
    df['sum']=df.apply(lambda x: x.sum(),axis=1)
    df.to_csv('new_sample.csv')
    for mode in range(0, 5, 1):
        np.savetxt('Spl_Para%i.txt'%mode, Para_mode[mode])
    pool.close()