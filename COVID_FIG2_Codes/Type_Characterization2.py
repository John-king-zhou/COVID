#provide functions for selecting the modes as defined
import numpy as np

def Peak(X,a,b):
    i=np.argmax(X)
    return max(X)>a and X[-1]<b

def PeakNum(X,threshold):
    count=0
    for i in range(1,len(X)-1,1):
        if (X[i]-X[i-1])>=0 and (X[i+1]-X[i])<=0 and X[i]>threshold:
            count+=1
    return count

def High(X,b):
    i=np.argmax(X)
    m=np.min(X[i:])
    return m>b

def Low(X,a):
    return max(X)<a

def within(x,a,b):
    return (min(x)>=a and max(x)<=b)

def WithinPhys(results):
    constraints={'v':1, 'Mphi':1, 'NK':1, 'Neut':1,
                 'T': 1, 'B':1, 'CD8>CD4': 1, 'ctk': 1}
    if not PeakNum(results[:,0],0.1)<3:
        constraints['v']=0
    if not within(results[:,4]+results[:,5],0,5):
        constraints['Mphi']=0
    if not within(results[:,6],0,0.5):
        constraints['NK']=0
    if not within(results[:,7],0,5):
        constraints['Neut']=0
    CD4T = np.sum(results[:, 8:17], axis=1)
    CD8T = np.sum(results[:, 17:21], axis=1)
    B = np.sum(results[:, 21:23], axis=1)
    if not within(CD4T+CD8T, 0, 8):
        constraints['T']=0
    if not within(B, 0, 5):
        constraints['B']=0
    # if np.max(CD8T)<np.max(CD4T):
    #     constraints['CD8>CD4']=0
    ylim=[500,200,50000,2000,1000,500]

    for j in range(24,29):
        if max(results[:,j])>ylim[j-24]:
            constraints['ctk'] = 0
    within_phys=(np.sum(list(constraints.values()))==8)
    return within_phys,constraints

def Type_Characterization(results):
    Vt=results[:,0]#virus time sequence
    IL6=results[:,26]
    if Vt[-1]<1e-6 and max(IL6)<1000:
        return 1
    if np.max(Vt)>1 and Vt[-1]<1e-6 and 1000<=max(IL6)<=2000:
        return 2
    if np.max(Vt)>1 and Vt[-1]<1e-6 and max(IL6)>2000:
        return 3
    if High(Vt,1) and max(IL6)>2000:
        return 4
    return 0