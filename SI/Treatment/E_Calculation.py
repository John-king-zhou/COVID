#provide functions for the calculation of gamma, e, e1, e2, e12
#args: results= 2D-array (len(time)*dimensions); Para=Parameters of the model, same as the input to ode integration
import numpy as np

def Hill(x,k,n):#Hill function
    return x**n/(k**n+x**n)

def load(results,Para):
    #number of cells of each kind in the lung
    global nCoV,If,H,D,APC_l,APC_u,NK,Neut,CD4Tn,CD4Ta,Th1,Th2,Th17,Tfh,Treg_a,Treg_r,CD4Tm,CD8Tn,CD8Ta
    global CTL,CD8Tm,Bgc,PB,Bm,IL2,IL4,IL6,IL10,TNF,IFN_g,Ig,A
    global k_APC_nCoV, k_APC_If, k_APC_rcr, k_NK_If, k_NK_APC, k_Neut_If, k_Neut_D, k_Neut_Th17, k_CD4_naive, k_CD4_mem, k_mem_CD4
    global k_Th1_CD4, k_Th2_CD4, k_Th17_CD4,k_Tfh_CD4, k_iTreg_CD4, k_nTreg_APC, k_CD8_naive, k_CD8_mem, k_CTL_CD8, k_mem_CD8
    global k_mem_CTL, k_GC_naive, k_GC_mem, k_PB, k_Bm
    global k_infect, d_v, k_c_1, k_c_2, k_c_3, k_c_4, k_k_1, k_k_2, k_k_3, k_k_4, k_k_5
    global h_APC_D, h_APC_TNF, h_APC_IFNg, h_NK_IL2,h_CD4_IL2, h_Th1_IFNg, h_Th2_IL4, h_Th17_IL6, h_Th17_Neut, h_Tfh_B
    global h_iTreg_IL10, h_CD8_IL2, h_CTL_Th1, h_CTL_IL2, h_CTL_IL6,h_Ig_IL4
    global r_H, r_APC, r_Treg_r, r_gc, m, d_APC_Treg, d_NK_Treg, d_CD4_Treg, d_CD8_Treg
    global d_If, d_H, d_D, d_APC_l, d_APC_u, d_NK, d_Neut, d_Th, d_Treg_a, d_Treg_r, d_CD4Tm, d_CTL, d_CD8Tm, d_GC, d_PB, d_Bm
    global p_IL2_0, p_IL2_1, p_IL2_2, p_IL2_3, p_IL2_4, p_IL4_0, p_IL4_1, p_IL6_0, p_IL6_1, p_IL6_2, p_IL6_3
    global p_IL10_0, p_IL10_1, p_IL10_2, p_TNF_0, p_TNF_1, p_TNF_2, p_TNF_3, p_IFN_g_0, p_IFN_g_1, p_IFN_g_2, p_IFN_g_3
    global p_Ig_1, p_Ig_2, c_IL2, c_IL4, c_IL6, c_IL10, c_TNF, c_IFN_g, c_Ig
    global K_IL2_1, K_IL2_2, K_IL2_3, K_IL2_4, K_IL2_5, K_IL2_6
    global K_IL4_1, K_IL4_2, K_IL4_3, K_IL6_1, K_IL6_2, K_IL6_3, K_IL10_1, K_IL10_2, K_IL10_3, K_IL10_4
    global K_TNF_1, K_IFNg_1, K_IFNg_2, K_IFNg_3, K_ACD4, K_ACD8, K_ACD8M, K_AB, K_m, K_nCoV_1, K_If_1, K_If_2, K_If_3
    global K_D_1, K_D_2, K_APC_1, K_APC_2, K_Neut_1, K_Th1_1, K_Th17_1, K_Tfh_1, K_GC_1
    global K_GC, APC0, NK0, Neut0, B0, N1, N_ex, t_CD4, t_CD8, g1, g2, g3, g4

    if len(results.shape)==1:
        nCoV=results[0]
        If=results[1]
        H=results[2]
        D=results[3]
        APC_l=results[4]
        APC_u=results[5]
        NK=results[6]
        Neut=results[7]
        CD4Tn=results[8]
        CD4Ta=results[9]
        Th1=results[10]
        Th2=results[11]
        Th17=results[12]
        Tfh=results[13]
        Treg_a=results[14]
        Treg_r=results[15]
        CD4Tm=results[16]
        CD8Tn=results[17]
        CD8Ta=results[18]
        CTL=results[19]
        CD8Tm=results[20]
        Bgc=results[21]
        PB=results[22]
        Bm=results[23]
        IL2=results[24]
        IL4=results[25]
        IL6=results[26]
        IL10=results[27]
        TNF=results[28]
        IFN_g=results[29]
        Ig=results[30]
        A=results[31]
    else:
        nCoV=results[:,0]
        If=results[:,1]
        H=results[:,2]
        D=results[:,3]
        APC_l=results[:,4]
        APC_u=results[:,5]
        NK=results[:,6]
        Neut=results[:,7]
        CD4Tn=results[:,8]
        CD4Ta=results[:,9]
        Th1=results[:,10]
        Th2=results[:,11]
        Th17=results[:,12]
        Tfh=results[:,13]
        Treg_a=results[:,14]
        Treg_r=results[:,15]
        CD4Tm=results[:,16]
        CD8Tn=results[:,17]
        CD8Ta=results[:,18]
        CTL=results[:,19]
        CD8Tm=results[:,20]
        Bgc=results[:,21]
        PB=results[:,22]
        Bm=results[:,23]
        IL2=results[:,24]
        IL4=results[:,25]
        IL6=results[:,26]
        IL10=results[:,27]
        TNF=results[:,28]
        IFN_g=results[:,29]
        Ig=results[:,30]
        A=results[:,31]

    #parameters
    k_a=Para[0:25]
    k_infect=Para[25]
    d_v=Para[26]
    k_c=Para[27:31]
    k_k=Para[31:36]
    h=Para[36:52]
    r=Para[52:56]
    m=Para[56]
    d_cell_Treg=Para[57:61]
    d_cell=Para[61:77]
    p_IL2=Para[77:82]
    p_IL4=Para[82:84]
    p_IL6=Para[84:88]
    p_IL10=Para[88:91]
    p_TNF=Para[91:95]
    p_IFN_g=Para[95:99]
    p_Ig=Para[99:101]
    c=Para[101:108]

    K_IL2=Para[108:114]
    K_IL4=Para[114:117]
    K_IL6=Para[117:120]
    K_IL10=Para[120:124]
    K_TNF_1=Para[124]
    K_IFN_g=Para[125:128]
    K_A=Para[128:132]
    K_m=Para[132]
    K_nCoV_1=Para[133]
    K_If=Para[134:137]
    K_D=Para[137:139]
    K_APC=Para[139:141]
    K_Neut_1=Para[141]
    K_Th1_1=Para[142]
    K_Th17_1=Para[143]
    K_Tfh_1=Para[144]
    K_GC_1=Para[145]
    K_GC=Para[146]
    APC0=Para[147]
    NK0=Para[148]
    Neut0=Para[149]
    B0=Para[150]
    N1=Para[151]
    N_ex=Para[152]
    t=Para[153:155]
    g=Para[155:159]

    k_APC_nCoV=k_a[0]
    k_APC_If=k_a[1]
    k_APC_rcr=k_a[2]
    k_NK_If=k_a[3]
    k_NK_APC=k_a[4]
    k_Neut_If=k_a[5]
    k_Neut_D=k_a[6]
    k_Neut_Th17=k_a[7]
    k_CD4_naive=k_a[8]
    k_CD4_mem=k_a[9]
    k_mem_CD4=k_a[10]
    k_Th1_CD4=k_a[11]
    k_Th2_CD4=k_a[12]
    k_Th17_CD4=k_a[13]
    k_Tfh_CD4=k_a[14]
    k_iTreg_CD4=k_a[15]
    k_nTreg_APC=k_a[16]
    k_CD8_naive=k_a[17]
    k_CD8_mem=k_a[18]
    k_CTL_CD8=k_a[19]
    k_mem_CD8=k_a[20]
    k_GC_naive=k_a[21]
    k_GC_mem=k_a[22]
    k_PB=k_a[23]
    k_Bm=k_a[24]

    k_c_1=k_c[0]
    k_c_2=k_c[1]
    k_c_3=k_c[2]
    k_c_4=k_c[3]

    k_k_1=k_k[0]
    k_k_2=k_k[1]
    k_k_3=k_k[2]
    k_k_4=k_k[3]
    k_k_5=k_k[4]

    h_APC_D=h[0]
    h_APC_TNF=h[1]
    h_APC_IFNg=h[2]
    h_NK_IL2=h[3]
    h_CD4_IL2=h[4]
    h_Th1_IFNg=h[5]
    h_Th2_IL4=h[6]
    h_Th17_IL6=h[7]
    h_Th17_Neut=h[8]
    h_Tfh_B=h[9]
    h_iTreg_IL10=h[10]
    h_CD8_IL2=h[11]
    h_CTL_Th1=h[12]
    h_CTL_IL2=h[13]
    h_CTL_IL6=h[14]
    h_Ig_IL4=h[15]

    r_H=r[0]
    r_APC=r[1]
    r_Treg_r=r[2]
    r_gc=r[3]

    d_APC_Treg=d_cell_Treg[0]
    d_NK_Treg=d_cell_Treg[1]
    d_CD4_Treg=d_cell_Treg[2]
    d_CD8_Treg=d_cell_Treg[3]

    d_If=d_cell[0]
    d_H=d_cell[1]
    d_D=d_cell[2]
    d_APC_l=d_cell[3]
    d_APC_u=d_cell[4]
    d_NK=d_cell[5]
    d_Neut=d_cell[6]
    d_Th=d_cell[7]
    d_Treg_a=d_cell[8]
    d_Treg_r=d_cell[9]
    d_CD4Tm=d_cell[10]
    d_CTL=d_cell[11]
    d_CD8Tm=d_cell[12]
    d_GC=d_cell[13]
    d_PB=d_cell[14]
    d_Bm=d_cell[15]

    p_IL2_0=p_IL2[0]
    p_IL2_1=p_IL2[1]
    p_IL2_2=p_IL2[2]
    p_IL2_3=p_IL2[3]
    p_IL2_4=p_IL2[4]

    p_IL4_0=p_IL4[0]
    p_IL4_1=p_IL4[1]

    p_IL6_0=p_IL6[0]
    p_IL6_1=p_IL6[1]
    p_IL6_2=p_IL6[2]
    p_IL6_3=p_IL6[3]

    p_IL10_0=p_IL10[0]
    p_IL10_1=p_IL10[1]
    p_IL10_2=p_IL10[2]

    p_TNF_0=p_TNF[0]
    p_TNF_1=p_TNF[1]
    p_TNF_2=p_TNF[2]
    p_TNF_3=p_TNF[3]

    p_IFN_g_0=p_IFN_g[0]
    p_IFN_g_1=p_IFN_g[1]
    p_IFN_g_2=p_IFN_g[2]
    p_IFN_g_3=p_IFN_g[3]

    p_Ig_1=p_Ig[0]
    p_Ig_2=p_Ig[1]

    c_IL2=c[0]
    c_IL4=c[1]
    c_IL6=c[2]
    c_IL10=c[3]
    c_TNF=c[4]
    c_IFN_g=c[5]
    c_Ig=c[6]

    K_IL2_1=K_IL2[0]
    K_IL2_2=K_IL2[1]
    K_IL2_3=K_IL2[2]
    K_IL2_4=K_IL2[3]
    K_IL2_5=K_IL2[4]
    K_IL2_6=K_IL2[5]

    K_IL4_1=K_IL4[0]
    K_IL4_2=K_IL4[1]
    K_IL4_3=K_IL4[2]

    K_IL6_1=K_IL6[0]
    K_IL6_2=K_IL6[1]
    K_IL6_3=K_IL6[2]

    K_IL10_1=K_IL10[0]
    K_IL10_2=K_IL10[1]
    K_IL10_3=K_IL10[2]
    K_IL10_4=K_IL10[3]

    K_IFNg_1=K_IFN_g[0]
    K_IFNg_2=K_IFN_g[1]
    K_IFNg_3=K_IFN_g[2]

    K_ACD4=K_A[0]
    K_ACD8=K_A[1]
    K_AB=K_A[2]
    K_mem=K_A[3]

    K_If_1=K_If[0]
    K_If_2=K_If[1]
    K_If_3=K_If[2]

    K_D_1=K_D[0]
    K_D_2=K_D[1]

    K_APC_1=K_APC[0]
    K_APC_2=K_APC[1]

    t_CD4=t[0]
    t_CD8=t[1]

    g1=g[0]
    g2=g[1]
    g3=g[2]
    g4=g[3]

def E(results,Para,period2=[0,0],period3=[0,0]):
    load(results,Para)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)
    IgA=np.multiply(Ig,A)
    e_clear=f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut+k_c_4*IgA
    e_kill=f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm+d_If
    if period2!=[0,0]:
        f_APC_anv[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        f_NK_eff[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        e_kill = f_APC_anv * k_k_1 * APC_l + f_APC_anv * k_k_2 * APC_u + f_NK_eff * k_k_3 * NK + k_k_4 * CTL + d_If
        e_clear = f_APC_anv * k_c_1 * APC_l + f_APC_anv * k_c_2 * APC_u + k_c_3 * Neut + k_c_4 * IgA
    if period3!=[0,0]:
        e_clear[int(period3[0]/0.1+1):min(int(period3[1]/0.1),results.shape[0]-1)]+=0.004*400
    e=e_clear*e_kill
    return e

def E_kill(results,Para,period2=[0,0],period3=[0,0]):
    load(results,Para)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)
    e_kill=f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_anv*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm+d_If
    if period2!=[0,0]:
        f_APC_anv[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        f_NK_eff[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        e_kill = f_APC_anv * k_k_1 * APC_l + f_APC_anv * k_k_2 * APC_u + f_NK_eff * k_k_3 * NK + k_k_4 * CTL + d_If
    return e_kill

def E_clear(results,Para,period2=[0,0],period3=[0,0]):
    load(results,Para)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    IgA=np.multiply(Ig,A)
    e_clear=f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut+k_c_4*IgA
    if period2!=[0,0]:
        f_APC_anv[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        f_NK_eff[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        e_clear = f_APC_anv * k_c_1 * APC_l + f_APC_anv * k_c_2 * APC_u + k_c_3 * Neut + k_c_4 * IgA
    if period3!=[0,0]:
        e_clear[int(period3[0]/0.1+1):min(int(period3[1]/0.1),results.shape[0]-1)]+=0.004*400
    return e_clear

def E1(results,Para,period2=[0,0],period3=[0,0]):#innate immunity
    load(results,Para)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)
    e_clear_1=f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut
    e_kill_1=f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+d_If
    if period2!=[0,0]:
        f_APC_anv[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        f_NK_eff[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        e_kill_1 = f_APC_anv * k_k_1 * APC_l + f_APC_anv * k_k_2 * APC_u + f_NK_eff * k_k_3 * NK + d_If
        e_clear_1 = f_APC_anv * k_c_1 * APC_l + f_APC_anv * k_c_2 * APC_u + k_c_3 * Neut
    e1=e_clear_1*e_kill_1
    return e1

def E2(results,Para,period2=[0,0],period3=[0,0]):#adaptive immunity
    load(results,Para)
    IgA=np.multiply(Ig,A)
    e_clear_2=k_c_4*IgA
    e_kill_2=k_k_4*CTL+k_k_5*CD8Tm+d_If
    if period3!=[0,0]:
        e_clear_2[int(period3[0]/0.1+1):min(int(period3[1]/0.1),results.shape[0]-1)]+=0.004*400
    e2=e_clear_2*e_kill_2
    return e2

def E12(results,Para,period2=[0,0],period3=[0,0]):
    load(results,Para)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)
    IgA=np.multiply(Ig,A)
    e_clear_1=f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut
    e_clear_2=k_c_4*IgA
    e_kill_1=f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK
    e_kill_2=+k_k_4*CTL+k_k_5*CD8Tm
    if period2!=[0,0]:
        f_APC_anv[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        f_NK_eff[int(period2[0]/0.1+1):min(int(period2[1]/0.1),results.shape[0]-1)]+=0.1
        e_kill_1 = f_APC_anv * k_k_1 * APC_l + f_APC_anv * k_k_2 * APC_u + f_NK_eff * k_k_3 * NK + d_If
        e_clear_1 = f_APC_anv * k_c_1 * APC_l + f_APC_anv * k_c_2 * APC_u + k_c_3 * Neut
    if period3!=[0,0]:
        e_clear_2[int(period3[0]/0.1+1):min(int(period3[1]/0.1),results.shape[0]-1)]+=0.004*400
    e_12=e_clear_1*e_kill_2+e_clear_2*e_kill_1
    return e_12
