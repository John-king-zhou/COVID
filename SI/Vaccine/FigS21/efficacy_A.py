import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import pandas as pd
import seaborn as sns
from Equation import func
from scipy.integrate import odeint

set2colors2 = ['#66c2a5', '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
flierprops = dict(marker='x', markerfacecolor='k', markersize=3,
                  linestyle='none', markeredgecolor='k')
ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
def get_results(Para,CD8Tm,Bm,A):
    time0=np.arange(0, 50, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159]/2, 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160]/2, 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    Ig0 = (1 + Para[51] * initial[25] / (initial[25] + Para[116])) * Bm * Para[100] / Para[107]
    initial[31] = A
    initial[30] = Ig0
    initial[23] = Bm
    initial[20] = CD8Tm
    initial[16] = 0.02
    results = odeint(func, initial, time0, args=(Para,))
    return results

def get_results_merge(args):
    return get_results(*args)
'''
eta-v0: simulation results
'''

r_seq = np.arange(0, 1, 0.2)
N8 = len(r_seq)
if __name__=='__main__':
    time = np.arange(0, 28, 0.1)
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores/3))
    A_seq=np.linspace(0,1,10)
    NA=len(A_seq)
    try:
        Protection=np.loadtxt('Protection_A.txt')
    except:
        sz=1000
        #Memory distribution
        mean = [-1, -1]
        std = [1, 1]
        Protection=np.zeros((NA,8))
        for mode in range(1,5,1):
            Parameter = np.loadtxt('GeoMean_para_Mode%i.txt' % mode)
            N=Parameter.shape[0]
            for j in range(NA):
                random_lgAb=np.random.normal(loc=mean[0],scale=std[0],size=sz)
                random_lgT=np.random.normal(loc=mean[1],scale=std[1],size=sz)
                random_points_Ab = 10**random_lgAb
                random_points_T = 10**random_lgT
                args=[]
                for i in range(sz):
                    args.append([Parameter,random_points_T[i],random_points_Ab[i],A_seq[j]])
                count=0
                Imm=0
                Mild=0
                print('A=%1.1f'%A_seq[j])
                for results in pool.imap(get_results_merge,args):
                        count+=1
                        ratio=count/sz
                        rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
                        rat_str=''.join(rat_str)
                        print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
                        v=results[:,0]
                        IL6=results[:,26]
                        if np.max(v) <= v[0]:
                            Imm+=1
                        if np.max(IL6)<2000:
                            Mild+=1
                Protection[j,mode-1]=Imm/sz
                Protection[j,mode+3]=Mild/sz
        Protection=np.array(Protection)
        np.savetxt('Protection_A.txt',Protection)
    pool.close()
    fig,ax=plt.subplots(1,1,figsize=(3.5,2))
    for mode in range(1,5,1):
        ax.plot(A_seq, Protection[:,mode-1]*100, c=ggcolors[mode], linewidth=1.5, ls='-')
        ax.plot(A_seq, Protection[:,mode+3]*100, c=ggcolors[mode], linewidth=1.5, ls='--')
    #ax.plot(logv_seq, np.mean(Protection[:, 0:4],axis=1)*100, c=ggcolors[0], linewidth=2.5)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,50,100])
    ax.set_ylim([-10,110])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Efficacy%',ha='center',va='center')
    ax.set_xlabel('Affinity',labelpad=2)
    labels=['Full Protection','Severe Prevention']
    ls=['-','--']
    l = [plt.plot([], [], color=ggcolors[0], ls=ls[i], label=labels[i])[0] for i in range(2)]
    ax.legend(handles=l, labels=labels, loc='upper right', bbox_to_anchor=([1.05, 1.35]), ncol=2, frameon=False,
               columnspacing=0.7, handlelength=1)
    fig.subplots_adjust(left=0.15,right=0.95,bottom=0.2,top=0.8,hspace=0.2)
    fig.savefig('efficacy_A.png')
    fig.savefig('efficacy_A.svg')
    plt.show()