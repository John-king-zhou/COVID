import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728']
set2colors = ['#66c2a5', '#fc8d62', '#e78ac3']
markers = ['o', '^', 's', 'd']
flierprops = dict(marker='x', markerfacecolor='k', markersize=3,
                  linestyle='none', markeredgecolor='k')
labels = ['Mode 1', '2', '3', '4', ]
vlabels = ['SARS-CoV-2', 'SARS', 'IAV']

df = pd.read_csv('Virus_Sample.csv')
v = np.array(df['vmax'])
df['logv']=np.log10(v)
# mode 1/2/3/4 ratio
fig1, ax1 = plt.subplots(1, 3, figsize=(8, 3))
for i in range(3):
    v=vlabels[i]
    df_v=df.loc[df.index[df['virus']==v]]
    proportion=[]
    for j in range(1,5,1):

        df_vj=df_v.loc[df_v.index[df_v['mode']==j]]
        #df_vj=df_v.iloc[ind]
        proportion.append(df_vj.shape[0])
    wedges, texts, autotexts=ax1[i].pie(proportion, normalize=True,
               colors=ggcolors,autopct='%1.1f%%', textprops={'color': 'w'})
    ax1[i].set_xlabel(vlabels[i])
fig1.legend(wedges, labels, framealpha=0, ncol=4, fontsize=12, bbox_to_anchor=(0.5,1), loc='upper center')
fig1.subplots_adjust(wspace=0,left=0,right=1)
fig1.savefig('virus_pie_chart.png',dpi=300)
fig1.savefig('virus_pie_chart.svg')

fig2,ax2=plt.subplots(1,4,figsize=(8,3.5))
columns=['logv','IL6max','Q','emax']
titles=['log$_{10}$V$_{max}$','IL-6$_{max}$','Q','$\epsilon_{max}$']
sFormatter1 = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
sFormatter1.set_powerlimits((0, 1))
for i in range(4):
    sns.boxplot(data=df,x='virus',y=columns[i],hue='mode',palette=ggcolors,ax=ax2[i],showfliers=False)
    ax2[i].set_xticklabels(vlabels,rotation=30,fontsize=10)
    ax2[i].set_ylabel('')
    ax2[i].set_xlabel('')
    ax2[i].set_title(titles[i])
    ax2[i].legend().remove()
ax2[3].yaxis.set_major_formatter(sFormatter1)
ax2[1].yaxis.set_major_formatter(sFormatter1)
fig2.subplots_adjust(left=0.05,right=0.95,bottom=0.18,top=0.9,wspace=0.25)
fig2.savefig('virus_boxplot_mode.png',dpi=300)

fig3,ax3=plt.subplots(1,4,figsize=(8,3.5))
for i in range(4):
    sns.boxplot(data=df,x='virus',y=columns[i],palette=set2colors,ax=ax3[i],showfliers=False)
    ax3[i].set_xticklabels(vlabels,rotation=30,fontsize=10)
    ax3[i].set_ylabel('')
    ax3[i].set_xlabel('')
    ax3[i].set_title(titles[i])
    ax3[i].legend().remove()
ax3[3].yaxis.set_major_formatter(sFormatter1)
ax3[1].yaxis.set_major_formatter(sFormatter1)
fig3.subplots_adjust(left=0.05,right=0.95,bottom=0.18,top=0.9,wspace=0.25)
fig3.savefig('virus_boxplot_together.png',dpi=300)
fig3.savefig('virus_boxplot_together.svg')

fig4,ax4=plt.subplots(1,3,figsize=(6,2.5))
Time = np.arange(0, 80, 0.1)
for i in range(3):
    v=vlabels[i]
    data=np.load('%s_ts.npy'%v)
    for j in range(3):
        y=data[j]
        if j==0:
            y[y<1e-4]=1e-4
            y=np.log10(y)
        ax4[j].plot(Time,np.mean(y,axis=0),color=set2colors[i],lw=2,label=v)
labels=['logV','IL-6 (pg/mL)','$\epsilon$']
for i in range(3):
    ax4[i].set_xlabel('time (days)')
    ax4[i].set_ylabel(labels[i])
    ax4[i].set_xlim([0,28])
    ax4[i].set_xticks([0,7,14,21,28])
    if i==1:
        ax4[i].yaxis.set_major_formatter(sFormatter1)
ax4[2].set_yticks([0,4,8])
fig4.subplots_adjust(bottom=0.18,top=0.8,left=0.1,right=0.95,wspace=0.3)
handles, labels = ax4[0].get_legend_handles_labels ()
fig4.legend(handles,labels,loc='upper center',bbox_to_anchor=(0.5,1),ncol=3,frameon=False)
ax4[1].set_yticks([0,1000,2000])
fig4.savefig('virus_comparison.png',dpi=300)
fig4.savefig('virus_comparison.svg')
plt.show()