'''
Use the data '3types_border_data.npy' generated from '_Border_3Types.py'.
figure 3D,E
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import scipy.interpolate as itp
def smooth_by_peak(x, y):
    sol_x = [x[0],]
    sol_y = [y[0],]
    for i in range(1,len(y)):
        if y[i]-y[i-1]!=0:
            sol_x.append((x[i]+x[i-1])/2)
            sol_y.append((y[i]+y[i-1])/2)
        else:
            sol_x.append(x[i])
            sol_y.append(y[i])

    return np.array(sol_x), np.array(sol_y)

def func_prot(x, g, cb, ct):
    '''
    y: lg(Tm)
    x: lg(Bm)
    k(Tm+ct)(Bm+cb)=gammg=3.6
    '''
    return np.log10(10**(g-np.log10(10**x+cb))-ct)


ggcolors = ['#808080', '#1F77B4', '#FF7F0E', '#D62728',
            '#2CA02C', '#4DBEEE', '#77AC30', '#9467BD']
set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
mode_name = ['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', ]
cm = LinearSegmentedColormap.from_list(
    'my_cm', [set2colors[3], 'w', set2colors[1]], N=100)
path = os.path.split(os.path.realpath(__file__))[0]

data = np.load('3types_border_dataA.npy', allow_pickle=True).item()#Ab only (exp)
data2 = np.load('3types_border_dataB.npy', allow_pickle=True).item()#Bm only (ss)
# --------- Figures ---------------------

CD8Tm_seq = [np.flip(data['CD8Tm']), data['CD8Tm'], data['CD8Tm']]
gamma = 3.6

# ------------- seperated borders for each modes in Bm,Tm space--------
figsp, axsp = plt.subplots(1, 4, figsize=(6, 2))
select_mode = [1, 2, 3, 4]
T8m_seq = CD8Tm_seq[0]
for im,mode in zip(range(4),select_mode):
    mode_list = np.array(data['mode'])
    indices = np.where(mode_list == mode)[0]
    Ab_data = np.array(data['bd_Bm'])
    Ab_data = Ab_data[indices]

    Ab_sum = [[] for i in range(len(T8m_seq))]
    for i in range(len(Ab_data)):
        Ab_tmp = Ab_data[i]
        for j in range(min([len(Ab_tmp), len(T8m_seq)])):
            Ab_sum[j].append(Ab_tmp[j])

    Ab_mean = np.array([np.average(np.array(seq)) for seq in Ab_sum])
    Ab_std = np.array([np.std(np.array(seq)) for seq in Ab_sum])
    #data2
    mode_list2 = np.array(data2['mode'])
    indices2 = np.where(mode_list2 == mode)[0]
    Ab_data2 = np.array(data2['bd_Bm'])
    Ab_data2 = Ab_data2[indices]

    Ab_sum2 = [[] for i in range(len(T8m_seq))]
    for i in range(len(Ab_data2)):
        Ab_tmp2 = Ab_data2[i]
        for j in range(min([len(Ab_tmp2), len(T8m_seq)])):
            Ab_sum2[j].append(Ab_tmp2[j])

    Ab_mean2 = np.array([np.average(np.array(seq)) for seq in Ab_sum2])
    Ab_std2 = np.array([np.std(np.array(seq)) for seq in Ab_sum2])

    axsp[im].plot(Ab_mean, T8m_seq, color=ggcolors[mode], lw=2)
    axsp[im].plot(Ab_mean2, T8m_seq, color=ggcolors[mode],lw=2, linestyle='--')
    axsp[im].fill_betweenx(T8m_seq, -Ab_std+Ab_mean, Ab_std+Ab_mean, color=ggcolors[mode], alpha=0.3,
                               edgecolor='w')

    axsp[im].set_xlim([0, 1300])
    axsp[im].set_xticks([0,500,1000])
    axsp[im].set_ylim([0, 1])
    axsp[im].set_yticks([0.5,1])
    axsp[im].set_yticklabels(['0.5','1'])
    axsp[im].text(x=1200,y=0.8,s=mode_name[mode-1],ha='right')
    if im in [0,]:
        axsp[im].set_yticklabels([0.5,1])
    else:
        axsp[im].set_yticklabels([])
axsp[0].plot([], [], color='k',lw=2, linestyle='-',label='exp')
axsp[0].plot([], [], color='k',lw=2, linestyle='--',label='ss')
figsp.legend(handlelength=2.5,frameon=False,ncol=2,bbox_to_anchor=(0.54,1),loc='upper center')
figsp.text(x=0.05,y=0.5,s='CD8+T$_M$ (t$^*$)',rotation=90,va='center',ha='center')
figsp.text(x=0.53,y=0.05,s='Ab (t$^*$)',rotation=0,va='center',ha='center')
figsp.subplots_adjust(left=0.12,right=0.95,bottom=0.23,top=0.8)
figsp.savefig('Modes_Border.svg')
figsp.savefig('Modes_Border.png')
fig2,ax2=plt.subplots(1,1,figsize=(2.2,2.2))
left, bottom, width, height = [0.6, 0.6, 0.25, 0.25]
ax3 = fig2.add_axes([left, bottom, width, height])
for dt, lstyle in zip([data, data2], ['-', '--']):
    log_ek_axis = np.linspace(np.min(np.log10(dt['bd_ek'])),
                              np.max(np.log10(dt['bd_ek'])),
                              25)

    for im, mode in zip(range(4), select_mode):
        log_ec_cmp = []
        mode_list = np.array(dt['mode'])
        indices = np.where(mode_list == mode)[0]
        ec_data = np.array(dt['bd_ec'])[indices]
        ek_data = np.array(dt['bd_ek'])[indices]

        for i in range(ec_data.shape[0]):
            x = np.log10(ec_data[i])
            y = np.log10(ek_data[i])
            y, x = smooth_by_peak(y, x)
            f = itp.interp1d(y, x, kind='linear',
                             bounds_error=False, fill_value=np.NaN)
            xxx = f(log_ek_axis)
            log_ec_cmp.append(xxx)

        log_ec_mean = np.nanmean(np.array(log_ec_cmp), axis=0)
        ax2.plot(log_ec_mean, log_ek_axis, color=ggcolors[mode], ls=lstyle)
        ax3.plot(log_ec_mean, log_ek_axis, color=ggcolors[mode], ls=lstyle)
x=np.arange(0.15,0.35,0.01)
import matplotlib.patches as patches
d=0.03
ax2.add_patch(
    patches.Rectangle((0.255, -0.045),d,d,facecolor='None',edgecolor='k',zorder=100))
ax2.plot(x,0.15-x,color='k')
ax2.text(x=0.13,y=-0.19,s='slope=-1',rotation=-45,fontsize=12)
ax3.plot([0.255,0.255+d],[-0.045+d,-0.045,],color='k')
ax3.text(x=0.262,y=-0.04,s='slope=-1',rotation=-45,ha='left',va='bottom',fontsize=9)
ax2.set_xlabel('log$_{10}\epsilon_c$',fontsize=13,labelpad=0)
ax2.set_ylabel('log$_{10}\epsilon_k$',fontsize=13,labelpad=0)
ax2.set_xticks([0.1,0.3,0.5])
ax2.set_yticks([-0.2,0,0.2])
ax2.set_aspect('equal')
ax2.set_xlim([0.1,0.5])
ax2.set_ylim([-0.2,0.2])
ax3.set_xlim([0.255,0.255+d])
ax3.set_ylim([-0.045,-0.045+d])
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_aspect('equal')
fig2.subplots_adjust(bottom=0.25,left=0.25)
fig2.savefig('ek-ec.svg')
fig2.savefig('ek-ec.png')
plt.show()
