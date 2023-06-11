import matplotlib.cm as cm
import matplotlib.mathtext as mathtext
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import interpolate


def Bm_defined_Ab(Para, Bm):
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159]/2, 0, 0, 0, 0, 0, 0,
               Para[54]/Para[70], 0, Para[160]/2, 0, 0, 0, 0, 0, 0, Para[77] /
               Para[101], Para[82] / Para[102],
               Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
               Para[104], Para[91]/Para[105], Para[95]/Para[106],
               0, 0]
    return (1 + Para[51] * initial[25] / (initial[25] + Para[116])) * Bm * Para[100] / Para[107]


# -----------General parameters-------------------------
set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
ggcolors = ['#808080', '#1F77B4', '#FF7F0E',
            '#D62728', '#4DBEEE', '#77AC30', '#9467BD']
mathtext.FontConstantsBase.sup1 = 0.4
time0 = np.arange(0, 50, 0.1)
ls = ['-', '--']
colorseq = [set2colors[1], set2colors[0], set2colors[3]]
path = os.path.split(os.path.abspath(__file__))[0]
cmap = cm.get_cmap('Blues_r')

# ---------------Select Mode---------------------------
mode = 3
Parameter = np.loadtxt(path+"\GeoMean_para_Mode%i.txt" % mode)
# -------------------------------------------------------------------------------
# ------------ Full Protection Surface ------------------------------------------
# -------------------------------------------------------------------------------
dt = np.load(path+'\\3D_AbBmTm_protect_m%d.npy' %
             mode, allow_pickle=True).item()
Bmseq, Abseq = dt['bm'], dt['Ab']
# ------------ 3D plot of Bm, Ab, Tm ----------------
XX, YY = np.meshgrid(Bmseq, Abseq)
ZZ_Tm = dt['Tm'][0]
for seq in dt['Tm'][1:]:
    ZZ_Tm = np.vstack((ZZ_Tm, seq))
ZZ_Tm = ZZ_Tm.T
# ZZ_Tm[ZZ_Tm<0]=0

# interpolate2d
x_seq = np.linspace(Bmseq[0], Bmseq[-1], 50)
y_seq = np.linspace(Abseq[0], Abseq[-1], 50)
XX_itp, YY_itp = np.meshgrid(x_seq,
                             y_seq)
f_zz = interpolate.interp2d(XX.flat, YY.flat, ZZ_Tm, kind='linear')
ZZ_plot = f_zz(x_seq, y_seq)

# figure settings
fig3d = plt.figure(figsize=(5, 4))
ax3d = fig3d.add_subplot(projection='3d')

# plot interpolated Z surface
ZZ_plot[Bm_defined_Ab(Parameter, XX_itp) >= YY_itp] = np.inf
# ax3d.plot_surface(XX_itp, YY_itp, ZZ_plot,
#                   facecolors=np.tile(np.array('green'),ZZ_plot.shape),
#                   alpha=0.2)

# plot original surface
ZZ_plot = np.copy(ZZ_Tm)
ZZ_plot[Bm_defined_Ab(Parameter, XX) >= YY] = np.NaN
ZZ_plot[ZZ_plot<0]=0
# ax3d.plot_wireframe(XX, YY, ZZ_plot, lw=1,color='gray'
#                     )
ax3d.plot_surface(XX, YY, ZZ_plot,
                  facecolors=np.tile(np.array('green'), ZZ_plot.shape),
                  alpha=0.2)

# -----------LIMIT Cases-------------
# Steady State ss
Bm_axis = np.arange(min(Bmseq), max(Bmseq), 0.01)
Ab_by_Bm = Bm_defined_Ab(Parameter, Bm_axis)
Bm_tmp = Bm_axis[Ab_by_Bm < max(Abseq)]
Ab_by_Bm = Ab_by_Bm[Ab_by_Bm < max(Abseq)]
Bm_tmp = Bm_tmp[Ab_by_Bm >= min(Abseq)]
Ab_by_Bm = Ab_by_Bm[Ab_by_Bm >= min(Abseq)]

zz_new = f_zz(Bm_tmp, Ab_by_Bm)
zz_line = np.array([zz_new[i, i] for i in range(len(zz_new))])

zz_line[zz_line<0]=0
ax3d.plot(Bm_tmp, Ab_by_Bm, zz_line,
          lw=1.5, alpha=0.8, ls='--',
          c='black', label='steady state (ss)')

# Exponentially decay
Bm_num_plot = [0, 4, 8]
lab_plot = ['exp decay',
            r'$B_M(t^*)$=%.1f' % Bmseq[Bm_num_plot[1]],
            r'$B_M(t^*)$=%.1f' % Bmseq[Bm_num_plot[2]]]
cmap_plot = [0, 0.3, 0.5]
mkr_plot = ['o', '^', 's']
for i in range(3):
    num_tmp = Bm_num_plot[i]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[:, num_tmp] for seq in [XX, YY, ZZ_Tm]]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[Bm_defined_Ab(
        Parameter, XX_tmp) <= YY_tmp] for seq in [XX_tmp, YY_tmp, ZZ_tmp]]
    ZZ_tmp[ZZ_tmp<0]=0
    ax3d.plot(XX_tmp, YY_tmp, ZZ_tmp,
              lw=1.5,
              label=lab_plot[i],
              color=cmap(cmap_plot[i]),
              marker=mkr_plot[i],
              zorder=100)

ax3d.set_xlim([0, 0.4])
ax3d.legend(loc='best')
ax3d.set_xlabel(r'$B_M(t^*)$')
ax3d.set_ylabel(r'$Ab(t^*)$')
ax3d.set_zlabel(r'$CD8+T_M(t^*)$')
ax3d.set_title('Full Proteciton Surface')
ax3d.view_init(23, 27)
fig3d.tight_layout()

# ----------Projection on Tm, Ab space-------------------
fig2d, ax2d = plt.subplots(figsize=(3.5, 3))
for i in range(3):
    num_tmp = Bm_num_plot[i]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[:, num_tmp] for seq in [XX, YY, ZZ_Tm]]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[Bm_defined_Ab(
        Parameter, XX_tmp) <= YY_tmp] for seq in [XX_tmp, YY_tmp, ZZ_tmp]]
    ZZ_tmp[ZZ_tmp<0]=0
    ax2d.plot(YY_tmp, ZZ_tmp,
              lw=1,
              label=lab_plot[i],
              color=cmap(cmap_plot[i]),
              marker=mkr_plot[i],
              markersize=4)
ax2d.plot(Ab_by_Bm, zz_line,
          lw=1.5, alpha=0.8, ls='--',
          c='black', label='steady state (ss)')

ax2d.legend()
ax2d.set_xlabel(r'$Ab(t^*)$')
ax2d.set_ylabel(r'$CD8+T_M(t^*)$')
ax2d.set_xlim([0,1000])
fig2d.tight_layout()

# -------------------------------------------------------------------------------
# ------------ Severe Prevention Surface ----------------------------------------
# -------------------------------------------------------------------------------
dt = np.load(path+'\\3D_AbBmTm_severePrevent_m%d.npy' %
             mode, allow_pickle=True).item()
Bmseq, Abseq = dt['bm'], dt['Ab']
# ------------ 3D plot of Bm, Ab, Tm ----------------
XX, YY = np.meshgrid(Bmseq, Abseq)
ZZ_Tm = dt['Tm'][0]
for seq in dt['Tm'][1:]:
    ZZ_Tm = np.vstack((ZZ_Tm, seq))
ZZ_Tm = ZZ_Tm.T
ZZ_Tm[ZZ_Tm < 0] = 0

# interpolate2d
x_seq = np.linspace(Bmseq[0], Bmseq[-1], 200)
y_seq = np.linspace(Abseq[0], Abseq[-1], 200)
XX_itp, YY_itp = np.meshgrid(x_seq,
                             y_seq)
f_zz = interpolate.interp2d(XX.flat, YY.flat, ZZ_Tm, kind='linear')
ZZ_plot = f_zz(x_seq, y_seq)

# figure settings
fig3d_s = plt.figure(figsize=(5, 4))
ax3d_s = fig3d_s.add_subplot(projection='3d')

# ZZ_plot[Bm_defined_Ab(Parameter, XX_itp)>=YY_itp]=np.inf
# ZZ_plot[ZZ_plot<0]=np.inf
# ax3d_s.plot_surface(XX_itp, YY_itp, ZZ_plot,
#                     cmap='Oranges_r',
#                     alpha=0.5)

# plot original surface
ZZ_plot = np.copy(ZZ_Tm)
ZZ_plot[Bm_defined_Ab(Parameter, XX) >= YY] = np.inf
# ZZ_plot[ZZ_plot<0]=np.NaN
ax3d_s.plot_surface(XX, YY, ZZ_plot,
                    facecolors=np.tile(np.array('purple'), ZZ_plot.shape),
                    alpha=0.2)

# Steady State ss
Bm_axis = np.arange(min(Bmseq), max(Bmseq)+0.005, 0.005)
Ab_by_Bm = Bm_defined_Ab(Parameter, Bm_axis)
Bm_tmp = Bm_axis[Ab_by_Bm < max(Abseq)]
Ab_by_Bm = Ab_by_Bm[Ab_by_Bm < max(Abseq)]
Bm_tmp = Bm_tmp[Ab_by_Bm >= min(Abseq)]
Ab_by_Bm = Ab_by_Bm[Ab_by_Bm >= min(Abseq)]

zz_new = f_zz(Bm_tmp, Ab_by_Bm)
zz_line = np.array([zz_new[i, i] for i in range(len(zz_new))])
Bm_tmp, Ab_by_Bm, zz_line = [seq[zz_line >= 0]
                             for seq in [Bm_tmp, Ab_by_Bm, zz_line]]

ax3d_s.plot(Bm_tmp, Ab_by_Bm, zz_line,
            lw=1.5, alpha=0.8, ls='--',
            c='black', label='steady state (ss)')

# Exponentially decay
Bm_num_plot = [0, 4, 8]
lab_plot = ['exp decay',
            r'$B_M(t^*)$=%.2f' % Bmseq[Bm_num_plot[1]],
            r'$B_M(t^*)$=%.2f' % Bmseq[Bm_num_plot[2]]]
cmap_plot = [0, 0.3, 0.5]
for i in range(3):
    num_tmp = Bm_num_plot[i]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[:, num_tmp] for seq in [XX, YY, ZZ_Tm]]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[Bm_defined_Ab(
        Parameter, XX_tmp) <= YY_tmp] for seq in [XX_tmp, YY_tmp, ZZ_tmp]]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[ZZ_tmp >= 0]
                              for seq in [XX_tmp, YY_tmp, ZZ_tmp]]
    ax3d_s.plot(XX_tmp, YY_tmp, ZZ_tmp,
                lw=1.5,
                label=lab_plot[i],
                color=cmap(cmap_plot[i]),
                marker=mkr_plot[i],
                fillstyle='none',
                zorder=100)

# figure settings
# ax3d_s.set_xlim(np.min(XX[ZZ_Tm>=0]),np.max(XX[ZZ_Tm>=0]))
# ax3d_s.set_ylim(np.min(YY[ZZ_Tm>=0]),np.max(YY[ZZ_Tm>=0]))
ax3d_s.set_xlim([0, 0.06])
ax3d_s.legend(loc='best')
ax3d_s.set_xlabel(r'$B_M(t^*)$')
ax3d_s.set_ylabel(r'$Ab(t^*)$')
ax3d_s.set_zlabel(r'$CD8+T_M(t^*)$')
ax3d_s.set_title('Severe Prevention Surface')
ax3d_s.view_init(23, 27)
fig3d_s.tight_layout()

# ----------Projection on Tm, Ab space-------------------
fig2d_s, ax2d_s = plt.subplots(figsize=(3.5, 3))
for i in range(3):
    num_tmp = Bm_num_plot[i]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[:, num_tmp] for seq in [XX, YY, ZZ_Tm]]
    XX_tmp, YY_tmp, ZZ_tmp = [seq[Bm_defined_Ab(
        Parameter, XX_tmp) <= YY_tmp] for seq in [XX_tmp, YY_tmp, ZZ_tmp]]
    ax2d_s.plot(YY_tmp, ZZ_tmp,
                lw=1,
                label=lab_plot[i],
                color=cmap(cmap_plot[i]),
                marker=mkr_plot[i],
                markersize=4,
                fillstyle='none')
ax2d_s.plot(Ab_by_Bm, zz_line,
            lw=1.5, alpha=0.8, ls='--',
            c='black', label='steady state (ss)')

ax2d_s.legend()
ax2d_s.set_xlabel(r'$Ab(t^*)$')
ax2d_s.set_ylabel(r'$CD8+T_M(t^*)$')
fig2d_s.tight_layout()


# ----------- Save figures ------------------
# fig3d.savefig('3d_full_protect.svg')
# fig2d.savefig('2d_full_protect.svg')
# fig3d_s.savefig('3d_severe_prevent.svg')
# fig2d_s.savefig('2d_severe_prevent.svg')

plt.show()
