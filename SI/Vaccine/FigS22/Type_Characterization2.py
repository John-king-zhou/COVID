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
    if not PeakNum(results[:,0],0.1)<3:
        return False
    if not within(results[:,4]+results[:,5],0,10):
        return False
    if not within(results[:,6],0,10):
        return False
    if not within(results[:,7],0,30):
        return False
    for i in range(8,24,1):
        if (not within(results[:,i],0,15)):
            return False
    CD4T = np.sum(results[:, 8:14], axis=1) + results[:, 16]
    CD8T = np.sum(results[:, 17:20], axis=1)
    if np.max(CD8T)<np.max(CD4T):
        return False
    for j in [24, 25, 27, 28, 29]:
        if max(results[:,j])>500:
            return False
    if not 5000>max(results[:,26]):
        return False
    # if not 3000>max(results[:,30]):
    #     return False
    return True

def Type_Characterization(results):
    Vt=results[:,0]#virus time sequence
    IL6=results[:,26]
    if np.max(Vt)>1 and Vt[-1]<1e-6 and 1000<=max(IL6)<=2000:
        return 2
    if np.max(Vt)>1 and Vt[-1]<1e-6 and max(IL6)>2000:
        return 3
    if High(Vt,1) and max(IL6)>2000:
        return 4
    if Vt[-1]<1e-6 and max(IL6)<1000:
        return 1
    return 0