import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
from Equation import func as func
from E_Calculation import *
import warnings

warnings.filterwarnings('error')

ggcolors=['#2CA02C','#1F77B4','#FF7F0E',]#'#D62728',]
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

def get_result(Para):
    time0=np.arange(0, 80, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    results = odeint(func, initial, time0, args=(Para,))
    CD4T = np.sum(results[:, 8:16], axis=1)
    CD8T = np.sum(results[:, 17:21], axis=1)
    CD4T[CD4T < 1e-4] = 1e-4
    CD8T[CD8T < 1e-4] = 1e-4
    Ntitre=NT(results,Para)
    results = np.vstack((results.T, CD4T, CD8T, Ntitre)).T
    return results

if __name__=='__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-1))
    #Ab comparison
    fig2,ax2=plt.subplots(2,1,figsize=(5,4))
    indices_Ab=np.arange(0,800,70)
    time_Ab=np.arange(0,80,7)
    df_Ab={'mode':[],'time':[],'Ab':[],'A':[],'NT':[]}
    try:
        df_Ab=pd.read_csv('verify_Ab.csv')
    except:
        for mode in [1,2,3]:
            print('mode=',mode)
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)#[0:20,:]
            data=[]
            N=Paras.shape[0]
            count=0
            for results in pool.imap(get_result,Paras):
                count+=1
                ratio=count/N
                rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
                rat_str=''.join(rat_str)
                print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
                data.append(results)
            data=np.array(data)
            Ab=data[:,indices_Ab,30]
            A=data[:,indices_Ab,31]
            Ntitre=data[:,indices_Ab,34]
            for i in range(len(time_Ab)):
                df_Ab['time']+=[int(time_Ab[i]) for j in range(Ab.shape[0])]
                df_Ab['mode']+=[int(mode) for j in range(Ab.shape[0])]
                df_Ab['Ab']+=list(Ab[:,i])
                df_Ab['A']+=list(A[:,i])
                df_Ab['NT']+=list(Ntitre[:,i])
        df_Ab=pd.DataFrame(df_Ab)

        df_Ab.to_csv('verify_Ab.csv')

    sns.barplot(data=df_Ab,x='time',y='Ab',hue='mode',palette=ggcolors,errwidth=1,capsize=0.2,ax=ax2[0])
    sns.barplot(data=df_Ab,x='time',y='NT',hue='mode',palette=ggcolors,errwidth=1,capsize=0.2,ax=ax2[1])
    ax2[1].get_legend().remove()
    ax2[0].set_xlabel('')
    ax2[0].set_ylabel('Ab titre')
    ax2[1].set_xlabel('time (days)')
    ax2[1].set_ylabel('Ab NS')
    fig2.subplots_adjust(left=0.2,bottom=0.15)

    fig2.subplots_adjust(bottom=0.15)

    fig2.savefig('Ab_verify.png',dpi=300)
    plt.show()