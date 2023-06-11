import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import pandas as pd
from Equation import func
from scipy.integrate import odeint
from Latin_Hypercube import LHSample

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
    try:
        Protection=np.loadtxt('Protection_g.txt')
    except:
        sz=100
        N=100
        #Memory distribution
        mean = [-1, -1]
        std = [1, 1]
        Protection=[]
        sample_indices4 = np.array([25, 61, 151, 52])
        for mode in range(1,5,1):
            Parameter = np.loadtxt('GeoMean_para_Mode%i.txt' % mode)
            LogParaBound = []
            for i in range(len(Parameter)):
                a = np.log10(Parameter[i])
                if i in sample_indices4:
                    LogParaBound.append([a - np.log10(5), a + np.log10(5)])
                else:
                    LogParaBound.append([a, a])
            LogPara_matrix = np.array(LHSample(len(LogParaBound), LogParaBound, N))
            Para_matrix = np.power(10, LogPara_matrix)
            print('Mode=%i'%mode)
            for j in range(N):
                ratio = j/N
                rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
                rat_str = ''.join(rat_str)
                print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
                Para=Para_matrix[j,:]
                gamma = Para[25] * Para[151] * Para[61] * Para[52] / Para[62]
                random_lgAb=np.random.normal(loc=mean[0],scale=std[0],size=sz)
                random_lgT=np.random.normal(loc=mean[1],scale=std[1],size=sz)
                random_points_Ab = 10**random_lgAb
                random_points_T = 10**random_lgT
                args=[]
                for i in range(sz):
                    args.append([Para,random_points_T[i],random_points_Ab[i],1])
                Imm=0
                Mild=0
                for results in pool.imap(get_results_merge,args):
                    v=results[:,0]
                    IL6=results[:,26]
                    if np.max(v) <= v[0]:
                        Imm+=1
                    if np.max(IL6)<2000:
                        Mild+=1
                Protection.append([gamma,Imm/sz,Mild/sz,mode])
        Protection=np.array(Protection)
        np.savetxt('Protection_g.txt',Protection)
    pool.close()
    fig,ax=plt.subplots(1,1,figsize=(3.5,2))
    Protection[:,1:3]*=100
    df=pd.DataFrame(Protection,columns=['gamma', 'Imm', 'Mild', 'Mode'])
    df['log_g']=np.log10(df.gamma)
    #sns.scatterplot(data=df,x='gamma',y='Imm',hue='Mode',palette=ggcolors[1:5],markers='o')
    for mode in range(1,5,1):
        inds=df.index[df.Mode==mode]
        df1=df.iloc[inds,:]
        rand_inds = np.random.choice(len(inds), 30)
        df2=df1.iloc[rand_inds,:]
        ax.scatter(df2.log_g,df2.Imm,edgecolor='w',lw=0.5,facecolor=ggcolors[mode],marker='o',alpha=0.8)
        ax.scatter(df2.log_g,df2.Mild,edgecolor=ggcolors[mode],facecolor='none',marker='^',alpha=0.8)
    ax.set_xlim([-0.2,2.1])
    ax.set_xticks([0,1,2])
    '''
    fitting full protection efficacy
    '''
    inds1 = df.index[df.Imm !=1]
    df1=df.iloc[inds1,:]
    z1 = np.polyfit(df1.log_g, df1.Imm, 1)
    p1 = np.poly1d(z1)
    print(p1)
    x2 = np.arange(-1, np.max(df1.log_g), 0.1)
    y2 = z1[0] * x2 + z1[1]
    y2[y2>100]=100
    ax.plot(x2, y2, color=ggcolors[0], ls='-', linewidth=2)
    '''
    fitting severe protection efficacy
    '''
    inds2 = df.index[df.Mild !=1]
    df2=df.iloc[inds2,:]
    z2 = np.polyfit(df2.log_g, df2.Mild, 1)
    p2 = np.poly1d(z2)
    print(p2)
    x3 = np.arange(-1, np.max(df2.log_g), 0.1)
    y3 = z2[0] * x3 + z2[1]
    y3[y3>100]=100
    ax.plot(x3, y3, color=ggcolors[0], ls='--', linewidth=2)
    # #ax.plot(logv_seq, np.mean(Protection[:, 0:4],axis=1)*100, c=ggcolors[0], linewidth=2.5)
    # ax.set_xticks([0,0.5,1])
    # ax.set_yticks([0,50,100])
    # ax.set_ylim([-10,110])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Efficacy%',ha='center',va='center')
    ax.set_xlabel('log$_{10}\gamma$',labelpad=0)
    labels=['Full Protection','Severe Prevention']
    markers=['o','^']
    filled=['k','w']
    l = [plt.scatter([], [], c='k', marker=markers[i], label=labels[i], facecolor=filled[i]) for i in range(2)]
    ax.legend(loc='upper center', bbox_to_anchor=([0.5, 1.35]), ncol=2, frameon=False,
               columnspacing=0.7, handlelength=1)
    fig.subplots_adjust(left=0.15,right=0.95,bottom=0.2,top=0.8,hspace=0.2)
    fig.savefig('efficacy_g.png')
    fig.savefig('efficacy_g.svg')
    plt.show()