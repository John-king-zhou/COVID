import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

patches=[mpatches.Patch(facecolor='w', edgecolor='k', linewidth=0.5, label='fixed',),
         mpatches.Patch(facecolor='w', edgecolor='k', linewidth=0.5, hatch='///', label='random')]
fixed=[np.loadtxt('Spl_Para%i.txt'%mode).shape[0] for mode in range(1,5,1)]
random=[np.loadtxt('rSpl_Para%i.txt'%mode).shape[0] for mode in range(1,5,1)]
ggcolors=['#2CA02C','#1F77B4','#FF7F0E', '#D62728',]
label=['fixed' for i in range(4)]+['random' for i in range(4)]
modes=[1,2,3,4,1,2,3,4]
df={'mode':modes,'size':fixed+random,'label':label}
fig,ax=plt.subplots(1,1,figsize=(4,4))
x=np.arange(1,5,1)
ax.bar(x-0.22,fixed,width=0.4,color=ggcolors,edgecolor='k',linewidth=1)
ax.bar(x+0.22,random,width=0.4,color=ggcolors,hatch='/',edgecolor='k',linewidth=1)
ax.set_xlabel('mode')
ax.set_ylabel('sample size (out of 2000)')
ax.set_yticks([0,100,200])
fig.legend(handles=patches, loc='center', ncol=2, bbox_to_anchor=(0.5, 0.94),
              columnspacing=1, frameon=False, fontsize=11)
fig.subplots_adjust(left=0.16)
fig.savefig('mode_distribution_comparison.png',dpi=300)
fig.savefig('mode_distribution_comparison.svg')
plt.show()
#sns.barplot(df,x='mode',y='size',hue='label',palette=ggcolors)