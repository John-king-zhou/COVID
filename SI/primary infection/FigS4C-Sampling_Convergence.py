#Sampling the parameters in the range given by logparabound.py, run GetBound.py before running this one
#used with multiprocessing for parallel computation
#output txt files as 2D-array (n*len(Para)) for further analysis
#Para0: all the parameter sets within physiological range; Para1~3: mode 1~3; Para4: asymptomatic patients
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Latin_Hypercube import LHSample
from Equation import func
from Type_Characterization2 import Type_Characterization,WithinPhys
import warnings
import time

warnings.filterwarnings('error')
ggcolors = ['#2CA02C', '#1F77B4', '#FF7F0E', '#D62728', ]
labels=['Mode 1','Mode 2','Mode 3','Mode 4',]
def get_mode(Para):
    Time=np.arange(0, 50, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    result=odeint(func, initial, Time, args=(Para,))
    if WithinPhys(result):
        mode=Type_Characterization(result)
        return mode,Para
    else:
        return 100,Para

if __name__=='__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores/3))
    N_list=[3*10,3*100,3*200,3*500,3*1000,3*2000,3*5000]
    try:
        All_Modes=np.load('Modes_Conv.npy')
    except:
        All_Modes=[]
        for trial in range(10):
            Modes=np.zeros((len(N_list),5))
            time0=time.perf_counter()
            LogParaBound=np.loadtxt('logparabound.txt')
            for i in range(len(N_list)):
                N=N_list[i]
                print('N=%i'%N)
                LogPara_matrix = np.array(LHSample(len(LogParaBound), LogParaBound, N))
                Para_matrix = np.power(10, LogPara_matrix)
                Para_mode = [[] for i in range(5)]
                count = 1
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
                    if count%(N//3)==0:
                        print('N=%i, processes completed:'%N,count,'Within_Phys_Range:',len(Para_mode[0]),
                              'Mode 1:%i/%i'%(len(Para_mode[1]),len(Para_mode[0])),
                              'Mode 2:%i/%i' % (len(Para_mode[2]), len(Para_mode[0])),
                              'Mode 3:%i/%i'%(len(Para_mode[3]),len(Para_mode[0])),
                              'Mode 4:%i/%i' % (len(Para_mode[4]), len(Para_mode[0])),
                              'time used:%1.1f s'%(time.perf_counter()-time0))
                M=len(Para_mode[1])+len(Para_mode[2])+len(Para_mode[3])+len(Para_mode[4])
                for mode in range(1,5,1):
                    Modes[i,mode-1]=len(Para_mode[mode])/M
                Modes[i,4]=M
            All_Modes.append(Modes)
        All_Modes=np.array(All_Modes)
        np.save('Modes_Conv.npy',All_Modes)
    pool.close()
    Mean=np.mean(All_Modes,axis=0)
    Std=np.std(All_Modes,axis=0)
    x=Mean[:,4]
    fig,ax=plt.subplots(1,1,figsize=(5,3))
    for i in [0,1,2,3]:
        ax.errorbar(x=N_list,y=Mean[:,i]*100,yerr=Std[:,i]*100,color=ggcolors[i],
                    ecolor=ggcolors[i],elinewidth=2,capsize=5,capthick=1.5,marker='o',label=labels[i])
    ax.set_xscale('log')
    ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1.03,1.18),ncol=4,columnspacing=1,handlelength=0.5)
    ax.set_xlabel('Sampling Size')
    ax.set_ylabel('Proportion %')
    fig.subplots_adjust(bottom=0.2)
    fig.savefig('Modes_Conv.svg')
    fig.savefig('Modes_Conv.png')
    plt.show()