 #Examine the effect of each drug (AntV, IFN-I, Ig and GC) used singly during differen stages
#early: PAD 0~7, middle: PAD 7~11, late: PAD 11+
#figure 4A
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from Equation import func_med
from scipy.integrate import odeint
import time

colors=plt.rcParams['axes.prop_cycle'].by_key()['color']

def Hill(x,k,n):
    if k==0:
        return 1
    else:
        return x**n/(k**n+x**n)

def get_results(Para0,period1,period2,period3,period4,mode=0):
    alpha=1
    epsilon_IFNI=1
    h_APC_IFNI=0
    h_NK_IFNI=0
    k_c_5=0.006
    Ig_ex=0
    beta_GC=1
    d_lymph_GC=0
    Para=np.hstack((Para0, alpha, epsilon_IFNI, h_APC_IFNI, h_NK_IFNI, k_c_5, Ig_ex, beta_GC, d_lymph_GC))
    funcc=func_med
    dt=0.1
    Time=np.arange(0,80,dt)

    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
               Para[91] / Para[105], Para[95] / Para[106], 0, 0]
    results0=odeint(funcc, initial, Time, args=(Para,))
    T=[0,7,14,80]
    times=[np.arange(T[i],T[i+1],dt) for i in range(3)]
    Traj=[]
    for i in range(3):
        if period1[0]<=T[i] and period1[1]>=T[i+1]:
            Para[161] = 0.8
        else:
            Para[161] = 1
        if period2[0]<=T[i] and period2[1]>=T[i+1]:
            Para[162] = 0.95
            Para[163] = 0.1
            Para[164] = 0.1
        else:
            Para[162] = 1
            Para[163] = 0
            Para[164] = 0
        if period3[0]<=T[i] and period3[1]>=T[i+1]:
            Para[166] = 400
        else:
            Para[166] = 0
        if period4[0]<=T[i] and period4[1]>=T[i+1]:
            Para[167] = 0.6
            Para[168] = 0.1
        else:
            Para[167]=1
            Para[168]=0
        if i==0:
            initial_i=initial*1
        else:
            initial_i=Traj[i-1][-1,:]
        results_i=odeint(funcc,initial_i,times[i],args=(Para,))
        Traj.append(results_i)
    results=np.vstack(Traj)
    v_final=results[-1,0]
    h_min=min(results[:,2])
    IL6_max=max(results[:,26])
    v_final_0=results0[-1,0]
    h_min_0=min(results0[:,2])
    IL6_max_0=max(results0[:,26])
    return v_final,h_min,IL6_max,v_final_0,h_min_0,IL6_max_0

def get_results_merge(args):
    return get_results(*args)

def score(vf,hm,IL6m):
    IL6c=2000
    hc=30
    vc=1
    S=(1+1*Hill(hm,hc,1))*(1+2*Hill(IL6c,IL6m,1))*(1+1*Hill(vc,vf,1))
    return S

def Update(Mat,i,Med,Sc):
    if Med==4:
        Mat[i,0]+=Sc
        Mat[i,1]+=Sc
    elif Med==5:
        Mat[i,1]+=Sc
        Mat[i,2]+=Sc
    elif Med==6:
        Mat[i,0]+=Sc
        Mat[i,1]+=Sc
        Mat[i,2]+=Sc
    else:
        Mat[i,Med-1]+=Sc

markers=['o','^','s','X']
ps=['early','middle','late']
vars=['e$_{max}$','nCoV$_{final}$','H$_{min}$','IL-6$_{max}$']
labels=['AntV','IFN-I','Ab','GC']

periods=[[0,0],[0,7],[7,14],[14,80],[0,14],[7,80],[0,80]]
if __name__=='__main__':
    fig,axes=plt.subplots(nrows=1,ncols=3,figsize=(10,4))
    ax=axes.flat
    fig2,axes2=plt.subplots(nrows=2,ncols=3,figsize=(8,6))
    ax2=axes2.flat
    cores=multiprocessing.cpu_count()
    pool=multiprocessing.Pool(processes=int(cores/3))
    t0=time.perf_counter()
    for mode in range(2, 5, 1):
        try:
            Sc_Mat=np.loadtxt('Sc_Mat%i.txt'%mode)
            Sc_Mat_std=np.loadtxt('Sc_Mat_std%i.txt'%mode)
        except:
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)
            Sc_Mat=np.zeros((4, 3))
            Sc_Mat_std=np.zeros((4,3))
            for i in range(0,4,1):#antv/ifn-I/Ig/gc
                for j in range(1,4,1):#early/middle/late
                    Sc=[]
                    T=[[0,0],[0,0],[0,0],[0,0]]
                    T[i]=periods[j]
                    args=[[Paras[k],T[0],T[1],T[2],T[3],mode] for k in range(Paras.shape[0])]
                    count=0
                    print('%s treatment duration: %s'%(labels[i],periods[j]))
                    for resultk in pool.imap(get_results_merge,args):
                        ratio = count / Paras.shape[0]
                        rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
                        rat_str = ''.join(rat_str)
                        print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
                        count+=1
                        v_final, h_min, IL6_max, v_final_0, h_min_0, IL6_max_0=resultk
                        delta_Q=(score(v_final, h_min, IL6_max)-score(v_final_0, h_min_0, IL6_max_0))/score(v_final_0, h_min_0, IL6_max_0)
                        Sc.append(delta_Q)
                    Sc=np.array(Sc)
                    Q1=np.percentile(Sc,25,interpolation='midpoint')
                    Q3=np.percentile(Sc,75,interpolation='midpoint')
                    IQR=Q3-Q1
                    indices=(Sc>Q3+1.5*IQR)+(Sc<Q1-1.5*IQR)
                    Sc2=np.mean(Sc[~indices])
                    Update(Sc_Mat,i,j,Sc2)
                    Sc[Sc>0.04]=0.04
                    Update(Sc_Mat_std,i,j,np.std(Sc[~indices]))
                print('Time used: %1.1f s'%(time.perf_counter()-t0),'mode %i'%mode,'drug',labels[i])
            np.savetxt('Sc_Mat%i.txt'%mode,Sc_Mat)
            np.savetxt('Sc_Mat_std%i.txt'%mode,Sc_Mat_std)
        Sc_Mat[Sc_Mat>0.04]=0.04
        Sc_Mat[Sc_Mat<-0.02]=-0.02
        for i in range(4):
            x=[1, 2, 3]
            y=Sc_Mat[i, :]
            yerr=Sc_Mat_std[i,:]
            ax[mode-2].plot(x, y, color=colors[i], linewidth=2)
            ax[mode-2].scatter(x, y, label=labels[i], color=colors[i], marker=markers[i], s=80)
            ax2[mode-2+3*(i//2)].errorbar(x=x, y=y, yerr=yerr, capsize=2, color=colors[i], linewidth=2)
            ax2[mode-2+3*(i//2)].scatter(x, y, label=labels[i], color=colors[i], marker=markers[i], s=80)
        ax[mode-2].set_xlim([0.8,3.2])
        ax[mode-2].hlines(0,0.8,3.2,zorder=-10,colors='k')
        ax[mode-2].set_xticks(x)
        ax[mode-2].set_xticklabels(ps)
        ax[mode-2].set_xlabel('Mode %i'%(mode),fontsize=14)
        ax[mode-2].tick_params(labelsize=13)
        ax2[mode-2].set_xlim([0.8,3.2])
        ax2[mode-2].hlines(0,0.8,3.2,zorder=-10,colors='k')
        ax2[mode-2].set_xticks(x)
        ax2[mode-2].set_xticklabels([])
        ax2[mode-2].tick_params(labelsize=13)
        ax2[mode+1].set_xlim([0.8,3.2])
        ax2[mode+1].hlines(0,0.8,3.2,zorder=-10,colors='k')
        ax2[mode+1].set_xticks(x)
        ax2[mode+1].set_xticklabels(ps)
        ax2[mode+1].set_xlabel('Mode %i'%mode,fontsize=14)
        ax2[mode+1].tick_params(labelsize=13)
    pool.close()
    for i in range(3):
        ax[i].set_yticks(np.arange(0,0.05,0.02))
        ax2[i].set_yticks(np.arange(0,0.05,0.02))
        ax2[i+3].set_yticks(np.arange(0,0.05,0.02))
        ax[i].set_ylim([-0.008,0.045])
        ax2[i].set_ylim([-0.008,0.045])
        ax2[i+3].set_ylim([-0.008,0.045])
        if i!=0:
            ax[i].set_yticklabels([])
            ax2[i].set_yticklabels([])
            ax2[i+3].set_yticklabels([])
        else:
            fig.text(x=0.085,y=0.87,s=r'10$^{-2}$',fontsize=12)
            ax[i].set_yticklabels([0,2,'> 4'])
            fig2.text(x=0.085,y=0.92,s=r'10$^{-2}$',fontsize=12)
            ax2[i].set_yticklabels([0,2,'> 4'])
            ax2[i+3].set_yticklabels([0,2,'> 4'])
    ax[0].set_ylabel(r'$\Delta$Q/Q',fontsize=13,labelpad=-2)
    ax[2].legend(ncol=4,loc='upper right',bbox_to_anchor=(1,1.15),framealpha=0,columnspacing=1,
                 handlelength=0.5,fontsize=11)
    fig.subplots_adjust(top=0.85,bottom=0.15,left=0.1,right=0.95,wspace=0.1)
    fig.savefig('Timing.svg')
    fig.savefig('Timing.png')
    handles, labels = ax[2].get_legend_handles_labels()
    ax2[0].set_ylabel(r'$\Delta$Q/Q',fontsize=13,labelpad=-2)
    ax2[3].set_ylabel(r'$\Delta$Q/Q',fontsize=13,labelpad=-2)
    ax2[2].legend(handles=handles,labels=labels,ncol=4,loc='upper right',bbox_to_anchor=(1,1.2),framealpha=0,columnspacing=1,
                 handlelength=0.5,fontsize=11)
    fig2.subplots_adjust(top=0.9,bottom=0.1,left=0.1,right=0.95,wspace=0.1,hspace=0.1)
    fig2.savefig('Timing_err.png')
    plt.show()