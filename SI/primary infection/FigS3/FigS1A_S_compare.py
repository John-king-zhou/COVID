# comparison of T supply between mode 1,2,3,4
import matplotlib.pyplot as plt
import numpy as np

# Parameters
ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#2CA02C','#4DBEEE','#77AC30','#9467BD']

# load mean trajectory
aver_traj = [[],[],[],[]]
for m in [1,2,3,4]:
    _, aver_traj[m-1] =  np.load('Data_M%s_S_comparison.npy' %m, allow_pickle=True)
time = np.arange(0, 50, 0.1)

# figure of T supply comparison
fig_s, ax_s = plt.subplots(1,2, figsize=(9,2.5))
ax_s = ax_s.flat
color_list=[ggcolors[i] for i in [4,1,2,3]]
for m in [1,2,3,4]:
    ax_s[0].plot(time, aver_traj[m-1][:len(time),-1], 
                 c = color_list[m-1], lw=2,
                 label = 'Mode %s' %m)
    ax_s[1].plot(time, aver_traj[m-1][:len(time),-2], 
                 c = color_list[m-1], lw=2,
                 label = 'Mode %s' %m)
ax_s[0].legend(loc='lower left', ncol=4, columnspacing=1,
               bbox_to_anchor = (-0.4,0.95))
ax_s[0].set_ylabel(r'$J_{S,CD8}$')
ax_s[1].set_ylabel(r'$J_{S,CD4}$')
ax_s[0].set_xlabel('Days')
ax_s[1].set_xlabel('Days')

ax_s[0].set_xticks([0,14,28,42,])
ax_s[1].set_xticks([0,14,28,42,])

ax_s[0].grid(alpha=0.5)
ax_s[1].grid(alpha=0.5)

plt.subplots_adjust(top=0.8,
                    bottom=0.2,
                    left=0.15,
                    right=0.95,
                    hspace=0.2,
                    wspace=0.4)
plt.savefig('S_comparison.png', dpi=200)
plt.savefig('S1A_S_comparison.svg', dpi=200)
plt.show()