#statistics of the fitting results, examine the relationship between maximum viral load,
#immune efficacy and recovery time
#figure 6, figure S17
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy.interpolate import interp1d
from scipy.stats import mannwhitneyu

warnings.filterwarnings('ignore')

set2colors = ['#66c2a5','#fc8d62','#a6d854','#e78ac3','#ffd92f','#e5c494','#b3b3b3']
ggcolors=['#1F77B4','#FF7F0E','#D62728','#9467BD','#77AC30',]
'''
Neant et al. PNAS 2021
'''
df_Neant=pd.read_csv('Neant_e.csv')
inds1=df_Neant.index[df_Neant.e_type==1].tolist()
df1_Neant=df_Neant.iloc[inds1,:]
inds2=df_Neant.index[df_Neant.e_type==2].tolist()
df2_Neant=df_Neant.iloc[inds2,:]
'''
Liu et al. Lancet ID 2020
'''
df_Liu=pd.read_csv('Liu_e.csv')
inds1=df_Liu.index[df_Liu.e_type==1].tolist()
df1_Liu=df_Liu.iloc[inds1,:]
inds2=df_Liu.index[df_Liu.e_type==2].tolist()
df2_Liu=df_Liu.iloc[inds2,:]
'''
Goyal et al. Science adv. 2020
'''
df_Goyal=pd.read_csv('Goyal_e.csv')
inds1=df_Goyal.index[df_Goyal.e_type==1].tolist()
df1_Goyal=df_Goyal.iloc[inds1,:]
inds2=df_Goyal.index[df_Goyal.e_type==2].tolist()
df2_Goyal=df_Goyal.iloc[inds2,:]

df11=pd.concat([df1_Neant,df1_Goyal,df1_Liu])
g = sns.JointGrid()
sns.scatterplot(data=df11,x='e',y='lgvmax',palette=sns.set_palette(ggcolors),ax=g.ax_joint,s=80)
sns.kdeplot(data=df11,x='e',palette=sns.set_palette(ggcolors),ax=g.ax_marg_x,
             linewidth=1,common_norm=False)
sns.kdeplot(data=df11,y='lgvmax',palette=sns.set_palette(ggcolors),ax=g.ax_marg_y,
             linewidth=1,common_norm=False)
ax=g.ax_joint
from scipy.stats import pearsonr
z1 = np.polyfit(df11.e,df11.lgvmax, 1)
p1 = np.poly1d(z1)
x2=np.arange(0,3.6,0.01)
y2=z1[0]*x2+z1[1]
r=pearsonr(df11.e,df11.lgvmax)
print(p1,'pearson=',r[0])
ax.plot(x2,y2,color='k',linewidth=2,zorder=10)
ax.text(x=2.5,y=6,s='p=%1.2f'%r[0])
ax.set_ylim([3,10.5])
ax.set_xlim([-0.2,3.7])
ax.set_xlabel(r'$\epsilon$')
ax.set_ylabel('log$_{10}$nCoV')
fig=plt.gcf()
fig.set_size_inches(3.5,3.5)
fig.subplots_adjust(bottom=0.15,top=0.95,left=0.15,right=0.95)
fig.savefig('early_e_vmax.png')
fig.savefig('early_e_vmax.svg')
fig2,ax2=plt.subplots(1,1)
bins=np.arange(0,3.5,0.5)
height=[]
std=[]
for i in range(len(bins)-1):
    indices=(df11.e>=bins[i]) & (df11.e<bins[i+1])
    height.append(np.nanmean(df11.loc[indices,'lgvmax']))
    std.append(np.nanstd(df11.loc[indices,'lgvmax']))
print(bins[:-1],height)
ax2.bar(x=bins[:-1],height=height,width=0.2,yerr=std,capsize=3)
ax2.set_xlabel(r'$\epsilon$')
ax2.set_ylabel('log$_{10}$nCoV')

df1_NCOV=pd.concat([df1_Neant,df1_Goyal,df1_Liu])
df2_NCOV=pd.concat([df2_Neant,df2_Goyal,df2_Liu])
print(df1_NCOV.shape,df2_NCOV.shape)
df1_NCOV['Virus']=['SARS-CoV-2' for i in range(df1_NCOV.shape[0])]
df2_NCOV['Virus']=['SARS-CoV-2' for i in range(df2_NCOV.shape[0])]

df1_NCOV=df1_NCOV.reset_index(drop=True)
th=1.1
inds11=df1_NCOV.index[df1_NCOV.e<th].tolist()
inds12=df1_NCOV.index[df1_NCOV.e>=th].tolist()
df11=df1_NCOV.iloc[inds11,:]
df11['low_e']=[1 for i in range(df11.shape[0])]
df12=df1_NCOV.iloc[inds12,:]
df12['low_e']=[0 for i in range(df12.shape[0])]
df1=pd.concat([df11,df12])

U1, p1 = mannwhitneyu(df11.lgvmax, df12.lgvmax)
print('e1<%1.1f v.s. e1>%1.1f: p='%(th,th), p1)


fig33,ax3=plt.subplots(1,1,figsize=(4,2))
fig3 = plt.figure(figsize=(6,3))
gs = fig3.add_gridspec(1, 9)
ax31=fig3.add_subplot(gs[:,1:4])
ax32=fig3.add_subplot(gs[:,5:])
sns.boxplot(data=df1,x='low_e',y='lgvmax',ax=ax3,palette=ggcolors[0:2],showfliers=False,boxprops={'facecolor':'None'})
sns.stripplot(data=df1,x='low_e',y='lgvmax',ax=ax3,dodge=True,jitter=0.2,palette=ggcolors[0:2])
ax3.set_xticks([0,1])
ax3.set_xticklabels(['$\epsilon\geq$%1.1f'%th,'$\epsilon<$%1.1f'%th])
ax3.xaxis.set_tick_params(pad=7)
ax3.set_xlabel('')
ax3.set_ylabel('log$_{10}$nCoV$_{max}$')
ax3.set_ylim([3,12.5])
ax3.set_yticks([4,7,10])
ax3.hlines(11,0,1,colors='k',lw=1)
if p1<0.001:
    s='*** p=%1.4f'%p1
elif p1<0.01:
    s='** p=%1.3f'%p1
elif p1<0.05:
    s='* p=%1.3f'%p1
else:
    s='p=%1.2f'%p1
ax3.text(x=0.5,y=11,s=s,fontsize=10,ha='center',va='bottom')
'''
declining phase, relationship between recovery time and e
'''
sns.scatterplot(data=df2_NCOV,x='e',y='recovery',color='k',ax=ax32)
df2_NCOV['loge']=np.log10(df2_NCOV.e)
df2_NCOV['log_rec']=np.log10(df2_NCOV.recovery)
df2_NCOV.replace(-np.inf,np.nan,inplace=True)
df2_NCOV_nan=df2_NCOV.dropna(axis=0,subset=['loge','log_rec'])
z1 = np.polyfit(df2_NCOV_nan.loc[:,'loge'],df2_NCOV_nan.loc[:,'log_rec'], 1)
p1 = np.poly1d(z1)
print(p1)
x2=np.arange(1.8,100,1)
y2=10**z1[1]*x2**z1[0]
ax32.plot(x2,y2,color='k',linewidth=1,zorder=-10)
ax32.text(60,40,s='T$\propto \epsilon^{%1.1f}$'%(z1[0]),ha='center',va='center',color='k',bbox = dict(boxstyle='round', facecolor='white', alpha=0.5))
ax32.set_xscale('log')
ax32.set_yscale('log')
ax32.set_xlim([1,500])
ax32.set_ylim([0.5,100])
ax32.set_xticks([1,10,100])
plt.xticks(fontsize=10)
ax32.set_yticks([1,10,50])
plt.yticks(fontsize=10)
ax32.set_xlabel('$\epsilon$',fontsize=13,labelpad=-2)
ax32.set_ylabel('Recovery Time')
fig3.subplots_adjust(bottom=0.2,left=0.1,right=0.95,wspace=0.3)
fig3.savefig('e-v.png')
fig3.savefig('e-v.svg')
fig33.subplots_adjust(bottom=0.2,top=0.9,left=0.2,right=0.95)
fig33.savefig('e11.svg')

age=['age < 65yrs',r'age$\geq$65yrs']

df1_Neant=df1_Neant.reset_index(drop=True)
inds21=df1_Neant.index[df1_Neant.old==age[0]].tolist()
inds22=df1_Neant.index[df1_Neant.old==age[1]].tolist()
df21=df1_Neant.iloc[inds21,:]
df22=df1_Neant.iloc[inds22,:]

U2, p2 = mannwhitneyu(df21.e, df22.e)
print('age<65yrs e1 v.s. age>=65yrs e1: p=', p2)

fig4,ax4=plt.subplots(1,1,figsize=(4,2))
sns.boxplot(data=df1_Neant,x='old',y='e',ax=ax4,palette=set2colors[0:2],showfliers=False,boxprops={'facecolor':'None'})
sns.stripplot(data=df1_Neant,x='old',y='e',ax=ax4,dodge=True,jitter=0.2,palette=set2colors[0:2])
# ax3.set_xticks([0,1])
# ax3.set_xticklabels(['$\epsilon<$1.5','$\epsilon\geq$1.5'])
ax4.set_xlabel('')
ax4.set_ylabel('Rising Phase $\epsilon$')
ax4.set_ylim([-0.5,6.2])
ax4.set_yticks([0,5])
ax4.hlines(5,0,1,colors='k',lw=1)
if p2<0.001:
    s='***'
elif p2<0.01:
    s='**'
elif p2<0.05:
    s='*'
else:
    s='p=%1.2f'%p2
ax4.text(x=0.5,y=5.6,s=s,fontsize=11,ha='center',va='center')
fig4.subplots_adjust(bottom=0.2,top=0.9,left=0.2,right=0.95)
fig4.savefig('Early_e_age.png')
fig4.savefig('Early_e_age.svg')
df2_Neant=df2_Neant.reset_index(drop=True)
inds31=df2_Neant.index[df2_Neant.old==age[0]].tolist()
inds32=df2_Neant.index[df2_Neant.old==age[1]].tolist()
df31=df2_Neant.iloc[inds31,:]
df32=df2_Neant.iloc[inds32,:]

U3, p3 = mannwhitneyu(df31.e, df32.e)
print('age<65yrs e2 v.s. age>=65yrs e2: p=', p3)

fig5,ax5=plt.subplots(1,1,figsize=(4,2))
sns.boxplot(data=df2_Neant,x='old',y='e',ax=ax5,palette=set2colors[0:2],showfliers=False,boxprops={'facecolor':'None'})
sns.stripplot(data=df2_Neant,x='old',y='e',ax=ax5,dodge=True,jitter=0.2,palette=set2colors[0:2])
# ax3.set_xticks([0,1])
# ax3.set_xticklabels(['$\epsilon<$1.5','$\epsilon\geq$1.5'])
ax5.set_xlabel('')
ax5.set_ylabel('Declining Phase $\epsilon$')
ax5.set_ylim([-0.5,37])
ax5.set_yticks([0,15,30])
ax5.hlines(30,0,1,colors='k',lw=1)
if p3<0.001:
    s='***'
elif p3<0.01:
    s='**'
elif p3<0.05:
    s='*'
else:
    s='p=%1.2f'%p3
ax5.text(x=0.5,y=33,s=s,fontsize=11,ha='center',va='center')
fig5.subplots_adjust(bottom=0.2,top=0.9,left=0.2,right=0.95)
fig5.savefig('Late_e_age.png')
fig5.savefig('Late_e_age.svg')
plt.show()