#data from Neat et al. PNAS 2021
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from FIT1 import FIT2, Imm1, Imm2, func1, func2
import warnings

warnings.filterwarnings('ignore')

set2colors = ['#66c2a5','#fc8d62','#a6d854','#e78ac3','#ffd92f','#e5c494','#b3b3b3']
age=['age < 65yrs',r'age$\geq$65yrs']
df=pd.read_csv('dataset_PNAS_SARSCOV2_655_anon.csv',header=1,sep=';',decimal=',')
inds=df.index[df.loc[:,'type']==1].tolist()
df=df.iloc[inds,:]
df=df.reset_index(drop=True)
'''
time
'''
df.loc[:,'time_monolix']=np.array(df.loc[:,'time_monolix'])-14
df.loc[:,'delai_monolix']=np.array(df.loc[:,'delai_monolix'])-14
indivs=list(set(df.loc[:,'ID']))
'''
preclude the ones with two or less data points
'''
new_indivs=[]
dfs=[]
for id in range(len(indivs)):
    indiv=indivs[id]
    inds=df.index[df.loc[:,'ID']==indiv].tolist()
    df_i=df.iloc[inds,:]
    time=np.array(df_i.loc[:,'time_monolix'])
    if len(time)<3:
        continue
    else:
        new_indivs.append(indiv)
        dfs.append(df_i)
print(len(dfs))
dfs=pd.concat(dfs)
dfs=dfs.reset_index(drop=True)
df2={'ID':[],'e':[],'e_type':[],'old':[],'lgvmax':[],'recovery':[]}
fig,axes=plt.subplots(11,12,figsize=(10,10))
axes=axes.flat
for i in range(len(new_indivs)):
    ratio = i / len(new_indivs)
    rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
    rat_str = ''.join(rat_str)
    print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
    id=new_indivs[i]
    ax=axes[i]
    inds=dfs.index[dfs.loc[:,'ID']==id].tolist()
    df_i=dfs.iloc[inds,:]
    old=df_i.loc[inds[0],'age_cat_cov']-1
    time=np.array(df_i.loc[:,'time_monolix'])
    t0=time[0]
    y = np.array(df_i.loc[:, 'y'])
    results=FIT2(time,y,id)
    e1,e2=np.nan,np.nan
    t_rec=0
    if results[0]:#when there is a rising phase
        popt=results[1]
        pcov=results[2]
        Imm3=results[3]
        time=time-time[0]
        e1=popt[1]*1
        e2=popt[2]*3
        time33 = np.arange(0, 50, 0.1)
        y33 = Imm3(time33, *popt)
        if len(time33[y33 == 2]) != 0:
            t_rec = time33[y33 == 2][0]-time33[np.argmax(y33)]
        else:
            t_rec = 50
        ax.plot(time33 + t0, y33, color=set2colors[old])
        ax.scatter(time+t0, y, s=30, marker='x', color=set2colors[old])
        ax.text(x=33,y=7.5,s='$\epsilon$1=%1.1f \n$\epsilon$2=%1.1f' % (e1, e2),fontsize=6,ha='center',va='center')
        ax.text(x=33, y=10, s='ID=%i' % (id),fontsize=6,ha='center',va='center')
    else:
        popt=results[2]
        pcov=results[3]
        if results[1] == 1:
            t_rec=0
            time = time - time[0]
            e2 = popt[1] * 1
            time11 = np.arange(0, 50, 0.1)
            y11=Imm1(time11, *popt)
            ax.plot(time11+t0, y11, color=set2colors[old])
        else:
            time = time - time[0]
            e2 = popt[1] * 3
            time22 = np.arange(0, 50, 0.1)
            y22=Imm2(time22, *popt)
            if len(time22[y22==2])!=0:
                t_rec=time22[y22==2][0]
            else:
                t_rec=50
            ax.plot(time22+t0, y22, color=set2colors[old])
        ax.scatter(time+t0, y, s=30, marker='x', color=set2colors[old])
        ax.text(x=33, y=7.5, s='$\epsilon$%i=%1.1f' % (results[1],e2),fontsize=6,ha='center',va='center')
        ax.text(x=33, y=10, s='ID=%i' % (id),fontsize=6,ha='center',va='center')
    ax.set_xlim([-2, 45])
    ax.set_xticks([0, 20, 40])
    ax.set_ylim([0.5, 13.5])
    ax.set_yticks([2,5,8,11])
    if i%12==0:
        ax.set_yticklabels([2,5,8,11])
    else:
        ax.set_yticklabels([])
    if i>=119:
        ax.set_xticklabels([0, 20, 40])
    else:
        ax.set_xticklabels([])
    if not np.isnan(e1):
        df2['ID']+=[id,]
        df2['e']+=[e1,]
        df2['e_type']+=[1,]
        df2['old']+=[age[old],]
        df2['lgvmax']+=[np.max(y),]
        df2['recovery']+=[t_rec,]
    if not np.isnan(e2):
        df2['ID']+=[id,]
        df2['e']+=[e2,]
        df2['e_type']+=[2,]
        df2['old']+=[age[old],]
        df2['lgvmax']+=[np.max(y),]
        df2['recovery']+=[t_rec,]
df2=pd.DataFrame(df2)
df2.to_csv('Neant_e.csv')
axes[-1].set_visible(False)
fig.subplots_adjust(bottom=0.1,top=0.98,left=0.1,right=0.98)
fig.text(0.5,0.05,s='Days post symptom onset',ha='center',va='center')
fig.text(0.05,0.5,s='Néant 2021 PNAS \n log$_{10}$nCoV',ha='center',va='center',rotation=90)
#fig.suptitle('Néant 2021 PNAS')
fig.savefig('Fit_All_Neant.png')
fig.savefig('Fit_All_Neant.svg')
plt.show()