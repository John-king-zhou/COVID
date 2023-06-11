#plotting the ensemble-averaged curves of Treg and Il-10 / TGF-beta
import numpy as np
import matplotlib.pyplot as plt


ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
indices=[14,15,27]

time = np.arange(0, 50, 0.1)
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(6,3))
axes = axes.flat
for mode in [1,2,3,4]:
    color=ggcolors[mode]
    Mean=np.loadtxt('Mean%i.txt'%mode)
    Std = np.loadtxt('Std%i.txt' % mode)
    print(Mean.shape)
    for i in range(3):
        mean=Mean[indices[i]]
        std=Std[indices[i]]

        ax=axes[i]
        ax.plot(time, mean, linewidth=3, c=color)
        ax.fill_between(time, mean-std, mean+std, facecolor=color, alpha=0.2)
labels=["Treg$^a$","Treg$^r$",r"IL-10/TGF-$\beta$"]
for i in range(3):
    ax=axes[i]
    ax.set_title(labels[i])
    ax.set_xlim([0,35])
    ax.set_xticks(np.arange(0,36,7))
axes[1].set_xlabel("time (days)")
fig.subplots_adjust(left=0.08,right=0.97,wspace=0.4,bottom=0.2)
fig.savefig('Immunosuppression.png',dpi=300)
plt.show()