#Class-based PCA method for dimension-reduction of clustered data
import numpy as np
import copy

#1st method of standardization (normally used)
#data= 2D array (# of samples * dimensions)
def standardize(data):
    meanvals=np.mean(data,axis=0)
    mean_removed=data-meanvals
    std_d=np.std(mean_removed,axis=0)
    return mean_removed/std_d

def standardize2(data):#an alternative method (standardize the data in the range [-1,1]
    meanvals=np.mean(data,axis=0)
    mean_removed=data-meanvals
    delta=np.max(mean_removed,axis=0)-np.min(mean_removed,axis=0)
    return mean_removed/delta

#data= 2D-array (number of samples * dimensions), n= int number of classes,
#indices= list indices of the start and end of each class,
#weight= array weight of each covariance matrix (length=n), d= int targeted dimensions
def CPCA(data,n,indices,weight,d,std,eigvalue=False):
    data2=copy.copy(data)
    if std==1:
        data2=standardize(data)
    if std==2:
        data2=standardize2(data)
    COM=[]
    CoVar=[]
    for i in range(n):
        classi=data2[indices[i]:indices[i+1],:]
        COM.append(np.mean(classi,axis=0))
        CoVar.append(np.cov(classi,rowvar=False))
    COM=np.array(COM)
    CoVar0=np.cov(COM,rowvar=False)
    Mat=CoVar0
    for i in range(n):
        Mat=Mat-weight[i]*CoVar[i]
    eigvals,eigvecs=np.linalg.eig(Mat)
    eigind=np.argsort(-eigvals)[0:d]
    eigu=eigvecs[:,eigind]
    eigvals=eigvals[eigind]
    data3=np.dot(data2,eigu)
    if eigvalue==True:
        return data3,eigvals,eigu
    else:
        return data3,eigu


