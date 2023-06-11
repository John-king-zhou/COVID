import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.integrate import odeint
from Equation import func as func
from E_Calculation import E_kill,E_clear
import scipy.interpolate as itp
import warnings

def get_result(Para):
    time0=np.arange(0, 80, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    results = odeint(func, initial, time0, args=(Para,))
    e_k = E_kill(results, Para)
    e_c = E_clear(results, Para)
    results = np.vstack((results.T, e_c, e_k)).T
    return results

def get_data(mode):
    T8m_seq = np.flip(data['CD8Tm'])
    mode_list = np.array(data['mode'])
    indices = np.where(mode_list == mode)[0]
    Ab_data = np.array(data['bd_Bm'])
    Ab_data = Ab_data[indices]

    Ab_sum = [[] for i in range(len(T8m_seq))]
    for i in range(len(Ab_data)):
        Ab_tmp = Ab_data[i]
        for j in range(min([len(Ab_tmp), len(T8m_seq)])):
            Ab_sum[j].append(Ab_tmp[j])
        # axsp[mode-1].plot(bm_tmp, T8m_seq, color=ggcolors[mode], alpha=0.1)

    Ab_mean = np.array([np.average(np.array(seq)) for seq in Ab_sum])
    Ab_std = np.array([np.std(np.array(seq)) for seq in Ab_sum])
    # data2
    mode_list2 = np.array(data2['mode'])
    indices2 = np.where(mode_list2 == mode)[0]
    Ab_data2 = np.array(data2['bd_Bm'])
    Ab_data2 = Ab_data2[indices]

    Ab_sum2 = [[] for i in range(len(T8m_seq))]
    for i in range(len(Ab_data2)):
        Ab_tmp2 = Ab_data2[i]
        for j in range(min([len(Ab_tmp2), len(T8m_seq)])):
            Ab_sum2[j].append(Ab_tmp2[j])
        # axsp[mode-1].plot(bm_tmp, T8m_seq, color=ggcolors[mode], alpha=0.1)

    Ab_mean2 = np.array([np.average(np.array(seq)) for seq in Ab_sum2])
    Ab_std2 = np.array([np.std(np.array(seq)) for seq in Ab_sum2])

    return T8m_seq,Ab_mean,Ab_mean2

warnings.filterwarnings('error')

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
set2colors2 = ['#66c2a5', '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3
labels=['Mode 1','Mode 2','Mode 3','Mode 4',]
time0=np.arange(0,80,0.1)
time=np.arange(0,28,0.1)
df={'ek':[],'ec':[],'Bm':[],'Ab':[],'CD8Tm':[],'v':[],'mode':[]}
list=[1,2,3,4]

for mode in list:
    print('mode=',mode)
    try:
        Mem=np.loadtxt('Mem%i.txt'%mode)
    except:
        Paras=np.loadtxt('Spl_Para%i.txt'%mode)
        N=Paras.shape[0]
        Mem=[]
        for k in range(Paras.shape[0]):
            ratio=k/N
            rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
            rat_str=''.join(rat_str)
            print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
            Para=Paras[k,:]
            results=get_result(Para)
            v=results[0:280,0]
            v[v<1e-4]=1e-4
            logv=np.mean(np.log10(v))
            final=results[-1,:]
            Mem.append([final[16],final[20],final[23],final[30],final[-2],final[-1],logv])
        Mem=np.array(Mem)
        np.savetxt('Mem%i.txt'%mode,Mem)
    df['CD8Tm']+=Mem[:,1].tolist()
    df['Bm']+=Mem[:,2].tolist()
    df['Ab']+=Mem[:,3].tolist()
    df['ec']+=Mem[:,4].tolist()
    df['ek']+=Mem[:,5].tolist()
    df['v']+=Mem[:,6].tolist()
    df['mode']+=[labels[mode-1] for i in range(Mem.shape[0])]
    #ax.scatter(Mem[:,3],Mem[:,4],color=ggcolors[mode],edgecolor='w')

data = np.load('3types_border_dataA.npy', allow_pickle=True).item()
data2 = np.load('3types_border_dataB.npy', allow_pickle=True).item()

df=pd.DataFrame(df)
colors=[ggcolors[i] for i in list]
fig,axes=plt.subplots(2,2,figsize=(5,5))
axes=axes.flat
count=0
df['protectedA']=False
df['protectedB']=False
for mode in [1,2,3,4]:
    x,y1,y2=get_data(mode)
    ax=axes[count]
    inds=df.index[(df['mode']==labels[mode-1])].tolist()
    df1=df.iloc[inds]
    f1=itp.interp1d(x=x,y=y1,kind='linear',fill_value='extrapolate')
    inds2=df1.index[f1(df1['CD8Tm'])<df1['Ab']]
    df.loc[inds2,'protectedA']=True
    f2=itp.interp1d(x=x,y=y2,kind='linear',fill_value='extrapolate')
    inds2=df1.index[f2(df1['CD8Tm'])<df1['Ab']]
    df.loc[inds2,'protectedB']=True

    sns.kdeplot(data=df1,x='Ab',y='CD8Tm',ax=ax,color=ggcolors[mode])
    ax.scatter(df1.Ab,df1.CD8Tm,facecolor=ggcolors[mode],edgecolor='k')
    #y1=f1(np.linspace(0,2,100))
    #x=np.linspace(0,2,100)
    ax.plot(y1,x,color='k',linewidth=1)
    #y2=f2(np.linspace(0,2,100))
    ax.plot(y2,x,color='k',linewidth=1,ls='--')
    ax.text(x=1900,y=1.4,s=labels[mode-1],ha='right',va='top',color=ggcolors[mode],fontsize=14)
    ax.set_ylim([0,1.5])
    ax.set_xlim([0,2000])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([0,1000,2000])
    ax.set_yticks([1,])
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    count+=1
ax.plot([],[],color='k',linewidth=1,label='exp')
#y2=f2(np.linspace(0,2,100))
ax.plot([],[],color='k',linewidth=1,ls='--',label='ss')
axes[2].set_xlabel('[Ab] ($\mu$g/mL)')
axes[2].set_ylabel('[CD8+T$_M$] (10$^6$ cells/mL)')
fig.legend(ncol=2,frameon=False,loc='upper center',bbox_to_anchor=(0.5,0.95))
fig.savefig('Immunogenicity.svg')

fig2,ax2=plt.subplots(1,1,figsize=(5,2))
m=np.zeros((4,4))
count=0
for mode in [1,2,3,4]:
    inds=df.index[df['mode'] == labels[mode-1]].tolist()
    df2=df.iloc[inds,:]
    i=0
    for p in [True,False]:
        indices=df2.index[df2['protectedA'] == p].tolist()
        m[i,count]=len(indices)
        i+=1
    for p in [True,False]:
        indices=df2.index[df2['protectedB'] == p].tolist()
        m[i,count]=len(indices)
        i+=1
    count+=1
y=np.zeros(4)
lgds=['protected','susceptible']
for i in range(2):
    x=m[i]/np.sum(m[0:2],axis=0)*100
    ax2.bar(np.arange(0,4,1)-0.15,x+y,color=set2colors2[i],zorder=-10*i,label=lgds[i],edgecolor='k',
            linestyle='-',width=0.3)
    y+=x
y=np.zeros(4)
for i in [2,3]:
    x=m[i]/np.sum(m[2:4],axis=0)*100
    ax2.bar(np.arange(0,4,1)+0.15,x+y,color=set2colors2[i-2],zorder=-10*i,edgecolor='k',
            linestyle='--',width=0.3)
    y+=x
ax2.set_ylabel('Protected %',labelpad=0)
xlabels2=['Mode %i'%(i+1) for i in range(4)]
ax2.set_xticks([0,1,2,3])
ax2.set_xticklabels(xlabels2)
ax2.legend(ncol=2,frameon=False,loc='upper right',bbox_to_anchor=(1.05,1.2))
fig2.subplots_adjust(top=0.87,left=0.15)
fig2.savefig('Protection_Ratio.svg')

fig3,ax3=plt.subplots(1,1)
sns.stripplot(data=df,x='mode',y='Ab')
plt.show()
