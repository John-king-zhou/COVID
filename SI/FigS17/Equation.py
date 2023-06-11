#ODE equations
import numpy as np

def Hill(x,k,n):
    return x**n/(k**n+x**n)

def logistic(x,K):
    return x*(1-x/K)

def func(x0,time,Para):
    #variables: number of cells of each kind in the lung/concentration of cytokines/Ig
    nCoV=x0[0]
    If=x0[1]
    H=x0[2]
    D=x0[3]
    APC_l=x0[4]
    APC_u=x0[5]
    NK=x0[6]
    Neut=x0[7]
    CD4Tn=x0[8]
    CD4Ta=x0[9]
    Th1=x0[10]
    Th2=x0[11]
    Th17=x0[12]
    Tfh=x0[13]
    Treg_a=x0[14]
    Treg_r=x0[15]
    CD4Tm=x0[16]
    CD8Tn=x0[17]
    CD8Ta=x0[18]
    CTL=x0[19]
    CD8Tm=x0[20]
    Bgc=x0[21]
    PB=x0[22]
    Bm=x0[23]
    IL2=x0[24]
    IL4=x0[25]
    IL6=x0[26]
    IL10=x0[27]
    TNF=x0[28]
    IFN_g=x0[29]
    Ig=x0[30]
    A=x0[31]

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

    f_APC_inf=1+h_APC_TNF*Hill(TNF,K_TNF_1,1)+h_APC_D*Hill(D,K_D_1,1)+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)

    A_CD4=f_APC_anv*Hill(APC_l,K_ACD4,2)
    A_CD8=f_APC_anv*Hill(APC_l,K_ACD8,2)
    A_B=f_APC_anv*Hill(APC_l,K_AB,2)

    k_Th1=k_Th1_CD4*(1+h_Th1_IFNg*Hill(IFN_g,K_IFNg_2,1))*Hill(K_IL4_1,IL4,1)\
          *Hill(K_IL10_1,IL10,1)

    k_Th2=k_Th2_CD4*(1+h_Th2_IL4*Hill(IL4,K_IL4_2,1))*Hill(K_IFNg_3,IFN_g,1)*Hill(K_IL10_1,IL10,1)

    k_Th17=k_Th17_CD4*(1+h_Th17_IL6*Hill(IL6,K_IL6_1,1)*Hill(IL10,K_IL10_2,1)+h_Th17_Neut*Hill(Neut,K_Neut_1,1))\
           *Hill(K_IL10_1,IL10,1)

    k_Tfh=k_Tfh_CD4*Hill(K_IL10_1,IL10,1)*(1+h_Tfh_B*Hill(Bgc,K_GC_1,1))

    k_iTreg=k_iTreg_CD4*(1+h_iTreg_IL10*Hill(IL10,K_IL10_3,1))*Hill(IL2,K_IL2_3,1)*Hill(K_IL6_2,IL6,1)

    k_CTL=k_CTL_CD8*(1+h_CTL_Th1*Hill(Th1,K_Th1_1,1)+h_CTL_IL2*Hill(IL2,K_IL2_6,1)
                                 +h_CTL_IL6*Hill(IL6,K_IL6_3,1))*Hill(K_IL10_4,IL10,1)

    dnCoVdt=N1*d_If*If-(f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut+k_c_4*A*Ig)*nCoV-d_v*Hill(nCoV,K_m,1)

    dIfdt=k_infect*nCoV*H-(f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm)*If-d_If*If

    dHdt=r_H-k_infect*nCoV*H-d_H*H

    dDdt=(f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm)*If+d_If*If-d_D*D

    # dAPCldt=(k_APC_nCoV*Hill(nCoV,K_nCoV_1,1)+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_Treg*Treg_a*APC_l-d_APC_l*APC_l
    #
    # dAPCudt=r_APC+k_APC_rcr*Hill(f_APC_eff*APC_l,K_APC_1,1)*APC0-\
    #         (k_APC_nCoV*Hill(nCoV,K_nCoV_1,1)+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_u*APC_u

    dAPCldt=(k_APC_nCoV*nCoV+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_Treg*Treg_a*APC_l-d_APC_l*APC_l

    dAPCudt=r_APC+k_APC_rcr*Hill(f_APC_inf*APC_l,K_APC_1,1)*APC0-\
            (k_APC_nCoV*nCoV+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_u*APC_u

    dNKdt=(k_NK_If*Hill(If,K_If_2,1)+k_NK_APC*Hill(f_APC_inf*APC_l,K_APC_2,1))*NK0-d_NK_Treg*Treg_a*NK-d_NK*NK

    dNeutdt=k_Neut_If*Hill(If,K_If_3,1)*Neut0+k_Neut_D*Hill(D,K_D_2,1)*Neut0+\
            k_Neut_Th17*Hill(Th17,K_Th17_1,1)*Neut0-d_Neut*Neut

    dCD4Tndt=-A_CD4*k_CD4_naive*CD4Tn

    dCD4Tadt=A_CD4*(1+h_CD4_IL2*Hill(IL2,K_IL2_2,1))*(k_CD4_naive*CD4Tn*2**g1/(t_CD4*g1)+k_CD4_mem*CD4Tm*2**g2/(t_CD4*g2))\
             -(k_Th1+k_Th2+k_Th17+k_Tfh+k_iTreg)*CD4Ta-k_mem_CD4*Hill(K_mem,APC_l,2)*CD4Ta-d_CD4_Treg*Treg_a*CD4Ta

    dTh1dt=k_Th1*CD4Ta-d_CD4_Treg*Treg_a*Th1-d_Th*Th1

    dTh2dt=k_Th2*CD4Ta-d_CD4_Treg*Treg_a*Th2-d_Th*Th2

    dTh17dt=k_Th17*CD4Ta-d_CD4_Treg*Treg_a*Th17-d_Th*Th17

    dTfhdt=k_Tfh*CD4Ta-d_CD4_Treg*Treg_a*Tfh-d_Th*Tfh

    dTreg_adt=k_iTreg*CD4Ta+k_nTreg_APC*A_CD4*Hill(IL2,K_IL2_4,1)*Treg_r-d_Treg_a*Treg_a

    dTreg_rdt=r_Treg_r-k_nTreg_APC*A_CD4*Hill(IL2,K_IL2_4,1)*Treg_r-d_Treg_r*Treg_r

    dCD4Tmdt=k_mem_CD4*Hill(K_mem,APC_l,2)*CD4Ta-k_CD4_mem*A_CD4*CD4Tm-d_CD4Tm*CD4Tm

    dCD8Tndt=-k_CD8_naive*A_CD8*CD8Tn

    dCD8Tadt=A_CD8*(k_CD8_naive*CD8Tn*(2**g3)/(t_CD8*g3)+k_CD8_mem*CD8Tm*(2**g4)/(t_CD8*g4))*(1+h_CD8_IL2*Hill(IL2,K_IL2_5,1))\
             -k_CTL*CD8Ta-d_CD8_Treg*Treg_a*CD8Ta-k_mem_CD8*Hill(K_mem,APC_l,2)*CD8Ta

    dCTLdt=k_CTL*CD8Ta-k_k_4*CTL*If/N_ex-d_CD8_Treg*Treg_a*CTL-d_CTL*CTL

    dCD8Tmdt=k_mem_CD8*Hill(K_mem,APC_l,2)*CD8Ta-k_CD8_mem*A_CD8*CD8Tm-d_CD8Tm*CD8Tm

    dBgcdt=A_B*(k_GC_naive*B0+k_GC_mem*Bm+r_gc*Hill(Tfh,K_Tfh_1,1)*Bgc)*(1-Bgc/K_GC)\
           -k_PB*A_B*Bgc-k_Bm*Bgc-d_GC*Bgc

    dPBdt=k_PB*A_B*Bgc-d_PB*PB

    dBmdt=k_Bm*Bgc-k_GC_mem*A_B*Bm-d_Bm*Bm

    dIL2dt=p_IL2_0+p_IL2_1*CD4Ta+p_IL2_2*CD8Ta+p_IL2_3*Th1+p_IL2_4*CTL-c_IL2*IL2

    dIL4dt=p_IL4_0+p_IL4_1*Th2-c_IL4*IL4

    dIL6dt=p_IL6_0+p_IL6_1*If+p_IL6_2*f_APC_inf*APC_l+p_IL6_3*Neut-c_IL6*IL6

    dIL10dt=p_IL10_0+p_IL10_1*Treg_a+p_IL10_2*Treg_r-c_IL10*IL10

    dTNFdt=p_TNF_0+p_TNF_1*If+p_TNF_2*f_APC_inf*APC_l+p_TNF_3*f_NK_eff*NK-c_TNF*TNF

    dIFNgdt=p_IFN_g_0+p_IFN_g_1*f_NK_eff*NK+p_IFN_g_2*Th1+p_IFN_g_3*CTL-c_IFN_g*IFN_g

    dIgdt=(1+h_Ig_IL4*Hill(IL4,K_IL4_3,1))*(p_Ig_1*PB+p_Ig_2*Bm)-c_Ig*Ig

    dAdt=m*Bgc*Hill(Tfh,K_Tfh_1,1)*(1-A)

    return np.array([dnCoVdt,dIfdt,dHdt,dDdt,dAPCldt,dAPCudt,dNKdt,dNeutdt,dCD4Tndt,dCD4Tadt,dTh1dt,dTh2dt,dTh17dt,
                     dTfhdt,dTreg_adt,dTreg_rdt,dCD4Tmdt,dCD8Tndt,dCD8Tadt,dCTLdt,dCD8Tmdt,dBgcdt,dPBdt,dBmdt,dIL2dt,
                     dIL4dt,dIL6dt,dIL10dt,dTNFdt,dIFNgdt,dIgdt,dAdt])

def func_IFNI(x0,time,Para):
    #variables: number of cells of each kind in the lung/concentration of cytokines/Ig
    nCoV=x0[0]
    If=x0[1]
    H=x0[2]
    D=x0[3]
    APC_l=x0[4]
    APC_u=x0[5]
    NK=x0[6]
    Neut=x0[7]
    CD4Tn=x0[8]
    CD4Ta=x0[9]
    Th1=x0[10]
    Th2=x0[11]
    Th17=x0[12]
    Tfh=x0[13]
    Treg_a=x0[14]
    Treg_r=x0[15]
    CD4Tm=x0[16]
    CD8Tn=x0[17]
    CD8Ta=x0[18]
    CTL=x0[19]
    CD8Tm=x0[20]
    Bgc=x0[21]
    PB=x0[22]
    Bm=x0[23]
    IL2=x0[24]
    IL4=x0[25]
    IL6=x0[26]
    IL10=x0[27]
    TNF=x0[28]
    IFN_g=x0[29]
    Ig=x0[30]
    A=x0[31]

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

    f_APC_inf=1+h_APC_TNF*Hill(TNF,K_TNF_1,1)+h_APC_D*Hill(D,K_D_1,1)
    f_APC_anv=1+h_APC_IFNg*Hill(IFN_g,K_IFNg_1,1)
    f_NK_eff=1+h_NK_IL2*Hill(IL2,K_IL2_1,1)

    A_CD4=f_APC_anv*Hill(APC_l,K_ACD4,2)
    A_CD8=f_APC_anv*Hill(APC_l,K_ACD8,2)
    A_B=f_APC_anv*Hill(APC_l,K_AB,2)

    k_Th1=k_Th1_CD4*(1+h_Th1_IFNg*Hill(IFN_g,K_IFNg_2,1))*Hill(K_IL4_1,IL4,1)\
          *Hill(K_IL10_1,IL10,1)

    k_Th2=k_Th2_CD4*(1+h_Th2_IL4*Hill(IL4,K_IL4_2,1))*Hill(K_IFNg_3,IFN_g,1)*Hill(K_IL10_1,IL10,1)

    k_Th17=k_Th17_CD4*(1+h_Th17_IL6*Hill(IL6,K_IL6_1,1)*Hill(IL10,K_IL10_2,1)+h_Th17_Neut*Hill(Neut,K_Neut_1,1))\
           *Hill(K_IL10_1,IL10,1)

    k_Tfh=k_Tfh_CD4*Hill(K_IL10_1,IL10,1)*(1+h_Tfh_B*Hill(Bgc,K_GC_1,1))

    k_iTreg=k_iTreg_CD4*(1+h_iTreg_IL10*Hill(IL10,K_IL10_3,1))*Hill(IL2,K_IL2_3,1)*Hill(K_IL6_2,IL6,1)

    k_CTL=k_CTL_CD8*(1+h_CTL_Th1*Hill(Th1,K_Th1_1,1)+h_CTL_IL2*Hill(IL2,K_IL2_6,1)
                                 +h_CTL_IL6*Hill(IL6,K_IL6_3,1))*Hill(K_IL10_4,IL10,1)

    dnCoVdt=N1*d_If*If-(f_APC_anv*k_c_1*APC_l+f_APC_anv*k_c_2*APC_u+k_c_3*Neut+k_c_4*A*Ig)*nCoV-d_v*Hill(nCoV,K_m,1)

    dIfdt=k_infect*nCoV*H-(f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm)*If-d_If*If

    dHdt=r_H-k_infect*nCoV*H-d_H*H

    dDdt=(f_APC_anv*k_k_1*APC_l+f_APC_anv*k_k_2*APC_u+f_NK_eff*k_k_3*NK+k_k_4*CTL+k_k_5*CD8Tm)*If+d_If*If-d_D*D

    # dAPCldt=(k_APC_nCoV*Hill(nCoV,K_nCoV_1,1)+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_Treg*Treg_a*APC_l-d_APC_l*APC_l
    #
    # dAPCudt=r_APC+k_APC_rcr*Hill(f_APC_eff*APC_l,K_APC_1,1)*APC0-\
    #         (k_APC_nCoV*Hill(nCoV,K_nCoV_1,1)+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_u*APC_u

    dAPCldt=(k_APC_nCoV*nCoV+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_Treg*Treg_a*APC_l-d_APC_l*APC_l

    dAPCudt=r_APC+k_APC_rcr*Hill(f_APC_inf*APC_l,K_APC_1,1)*APC0-\
            (k_APC_nCoV*nCoV+k_APC_If*Hill(If,K_If_1,1))*APC_u-d_APC_u*APC_u

    dNKdt=(k_NK_If*Hill(If,K_If_2,1)+k_NK_APC*Hill(f_APC_inf*APC_l,K_APC_2,1))*NK0-d_NK_Treg*Treg_a*NK-d_NK*NK

    dNeutdt=k_Neut_If*Hill(If,K_If_3,1)*Neut0+k_Neut_D*Hill(D,K_D_2,1)*Neut0+\
            k_Neut_Th17*Hill(Th17,K_Th17_1,1)*Neut0-d_Neut*Neut

    dCD4Tndt=-A_CD4*k_CD4_naive*CD4Tn

    dCD4Tadt=A_CD4*(1+h_CD4_IL2*Hill(IL2,K_IL2_2,1))*(k_CD4_naive*CD4Tn*2**g1/(t_CD4*g1)+k_CD4_mem*CD4Tm*2**g2/(t_CD4*g2))\
             -(k_Th1+k_Th2+k_Th17+k_Tfh+k_iTreg)*CD4Ta-k_mem_CD4*Hill(K_mem,APC_l,2)*CD4Ta-d_CD4_Treg*Treg_a*CD4Ta

    dTh1dt=k_Th1*CD4Ta-d_CD4_Treg*Treg_a*Th1-d_Th*Th1

    dTh2dt=k_Th2*CD4Ta-d_CD4_Treg*Treg_a*Th2-d_Th*Th2

    dTh17dt=k_Th17*CD4Ta-d_CD4_Treg*Treg_a*Th17-d_Th*Th17

    dTfhdt=k_Tfh*CD4Ta-d_CD4_Treg*Treg_a*Tfh-d_Th*Tfh

    dTreg_adt=k_iTreg*CD4Ta+k_nTreg_APC*A_CD4*Hill(IL2,K_IL2_4,1)*Treg_r-d_Treg_a*Treg_a

    dTreg_rdt=r_Treg_r-k_nTreg_APC*A_CD4*Hill(IL2,K_IL2_4,1)*Treg_r-d_Treg_r*Treg_r

    dCD4Tmdt=k_mem_CD4*Hill(K_mem,APC_l,2)*CD4Ta-k_CD4_mem*A_CD4*CD4Tm-d_CD4Tm*CD4Tm

    dCD8Tndt=-k_CD8_naive*A_CD8*CD8Tn

    dCD8Tadt=A_CD8*(k_CD8_naive*CD8Tn*(2**g3)/(t_CD8*g3)+k_CD8_mem*CD8Tm*(2**g4)/(t_CD8*g4))*(1+h_CD8_IL2*Hill(IL2,K_IL2_5,1))\
             -k_CTL*CD8Ta-d_CD8_Treg*Treg_a*CD8Ta-k_mem_CD8*Hill(K_mem,APC_l,2)*CD8Ta

    dCTLdt=k_CTL*CD8Ta-k_k_4*CTL*If/N_ex-d_CD8_Treg*Treg_a*CTL-d_CTL*CTL

    dCD8Tmdt=k_mem_CD8*Hill(K_mem,APC_l,2)*CD8Ta-k_CD8_mem*A_CD8*CD8Tm-d_CD8Tm*CD8Tm

    dBgcdt=A_B*(k_GC_naive*B0+k_GC_mem*Bm+r_gc*Hill(Tfh,K_Tfh_1,1)*Bgc)*(1-Bgc/K_GC)\
           -k_PB*A_B*Bgc-k_Bm*Bgc-d_GC*Bgc

    dPBdt=k_PB*A_B*Bgc-d_PB*PB

    dBmdt=k_Bm*Bgc-k_GC_mem*A_B*Bm-d_Bm*Bm

    dIL2dt=p_IL2_0+p_IL2_1*CD4Ta+p_IL2_2*CD8Ta+p_IL2_3*Th1+p_IL2_4*CTL-c_IL2*IL2

    dIL4dt=p_IL4_0+p_IL4_1*Th2-c_IL4*IL4

    dIL6dt=p_IL6_0+p_IL6_1*If+p_IL6_2*f_APC_inf*APC_l+p_IL6_3*Neut-c_IL6*IL6

    dIL10dt=p_IL10_0+p_IL10_1*Treg_a+p_IL10_2*Treg_r-c_IL10*IL10

    dTNFdt=p_TNF_0+p_TNF_1*If+p_TNF_2*f_APC_inf*APC_l+p_TNF_3*f_NK_eff*NK-c_TNF*TNF

    dIFNgdt=p_IFN_g_0+p_IFN_g_1*f_NK_eff*NK+p_IFN_g_2*Th1+p_IFN_g_3*CTL-c_IFN_g*IFN_g

    dIgdt=(1+h_Ig_IL4*Hill(IL4,K_IL4_3,1))*(p_Ig_1*PB+p_Ig_2*Bm)-c_Ig*Ig

    dAdt=m*Bgc*Hill(Tfh,K_Tfh_1,1)*(1-A)

    return np.array([dnCoVdt,dIfdt,dHdt,dDdt,dAPCldt,dAPCudt,dNKdt,dNeutdt,dCD4Tndt,dCD4Tadt,dTh1dt,dTh2dt,dTh17dt,
                     dTfhdt,dTreg_adt,dTreg_rdt,dCD4Tmdt,dCD8Tndt,dCD8Tadt,dCTLdt,dCD8Tmdt,dBgcdt,dPBdt,dBmdt,dIL2dt,
                     dIL4dt,dIL6dt,dIL10dt,dTNFdt,dIFNgdt,dIgdt,dAdt])