import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

original=[np.loadtxt('Spl_Para%i.txt'%mode).shape[0] for mode in range(1,5,1)]
test1=[np.loadtxt('./Test_Criteria/Spl_Para%i_%i.txt' % (0,mode)).shape[0] for mode in range(1,5,1)]
test2=[np.loadtxt('./Test_Criteria/Spl_Para%i_%i.txt' % (1,mode)).shape[0] for mode in range(1,5,1)]
ggcolors=['#2CA02C','#1F77B4','#FF7F0E', '#D62728',]

fig,ax=plt.subplots(1,1,figsize=(4,3))
x=np.linspace(-0.3,0.3,4)
ax.bar(x,original,width=0.2,color=ggcolors,edgecolor='k',linewidth=1)
ax.bar(x+1,test1,width=0.2,color=ggcolors,edgecolor='k',linewidth=1)
ax.bar(x+2,test2,width=0.2,color=ggcolors,edgecolor='k',linewidth=1)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["original","a=800\nb=1600","a=1200\nb=2400",])
ax.set_ylabel('sample size (out of 2000)')
ax.set_yticks([0,100,200,300])
fig.subplots_adjust(left=0.16,bottom=0.24)
fig.savefig('./Test_Criteria/test_criteria_distrib.png',dpi=300)

plt.show()
#sns.barplot(df,x='mode',y='size',hue='label',palette=ggcolors)