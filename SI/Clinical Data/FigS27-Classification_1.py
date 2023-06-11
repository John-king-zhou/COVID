#categorize patients' data into class 1, 2 and 3
#generate Classification1.xls file and figure 3A
#class 1 = max(IL-6)<170 & mean(e*)>0.13
#class 2 = max(IL-6)>170 & mean(e*)>0.13
#class 3 = max(IL-6)>170 & mean(e*)<0.13
import numpy as np
import xlwt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

workbook = xlwt.Workbook()
worksheet = workbook.add_sheet('Classification')
worksheet.write(0,0,'Number')
worksheet.write(0,1,'Class')
worksheet.write(0,2,'Severity')

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
def get_type(v,w):#v=IL6, w=e
    t1=170
    t2=0.15
    IL6max=np.nanmax(v)
    mean_e=np.nanmean(w)
    #print(mean_e)
    #mean_w=np.trapz(y=wt,x=t)/(t[-1]-t[0])
    if IL6max<100 and mean_e>t2:
        return 1
    elif IL6max>100 and IL6max<t1 and mean_e>t2:
        return 2
    elif IL6max>t1 and mean_e>t2:
        return 3
    elif IL6max>t1 and mean_e<t2:
        return 4
    else:
        return 0

Severity=np.zeros(3)
Class=np.zeros(5)
indivs=[]
data=[[] for i in range(3)]
CStat=['Mild/Moderate','Severe','Critical']
CLabel=['Mild/\nModerate','Severe','Critical']

row_n=0

excl = pd.ExcelFile('Time_Course.xlsx')
excl2 = pd.read_excel(excl,sheet_name=None,index_col=None)
keys = list(excl2.keys())
for filename in keys:
    s=0
    c=0
    df = pd.read_excel(excl, sheet_name=filename,header=0,index_col=None)
    data = df.iloc[1:,3:]
    #data=np.array(data).T
    row_n+=1
    s=int(df.loc[0,'Clinical Condition'])
    if s!=0:
        Severity[s-1]+=1
    IL6=np.array(data['IL-6'])
    E=np.array(data['$\epsilon$*'])
    c=get_type(IL6,E)
    Class[c]+=1
    indivs.append([c,s])
    worksheet.write(row_n, 0, filename)
    worksheet.write(row_n, 1, '%i'%c)
    worksheet.write(row_n, 2, CStat[s-1])
    worksheet.write(row_n, 3, np.mean(E[~np.isnan(E)]))
    worksheet.write(row_n, 4, np.max(IL6[~np.isnan(IL6)]))
workbook.save('Classification1.xls')
fig2=plt.figure(figsize=(5,5))
ax2=fig2.add_subplot(121)
sorted_indivs=sorted(indivs,key=lambda x: (x[0], x[1]))
sorted_indivs=np.array(sorted_indivs)
x=sorted_indivs[:,1]
for j in range(0,5):
    x=sorted_indivs[sorted_indivs[:,0]==j,1]
    print('Class %i'%j,'Mild:',sum(x==1),'Severe:',sum(x==2),'Critical:',sum(x==3))
sorted_indivs=np.vstack((sorted_indivs[:,0],sorted_indivs.T)).T
cmaps = colors.LinearSegmentedColormap.from_list('mylist',ggcolors,N=8)
sorted_indivs[:,2]+=4
print(sorted_indivs)
g1=ax2.imshow(sorted_indivs,cmap=cmaps,origin='lower',aspect='auto')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('Modes',size=13)
Types=['Others','Mode 1','Mode 2','Mode 3','Mode 4']
b=0
for i in range(0,5):
    if i==0:
        ax2.text(s=Types[i], x=0.5, y=b + Class[i] / 2 + 0, horizontalalignment='center',
                 verticalalignment='center', c='white', size=12, weight='bold')
    else:
        ax2.text(s=Types[i], x=0.5,y=b+Class[i]/2+0.9, horizontalalignment='center',
                 verticalalignment='center', c='white',size=12,weight='bold')
    if Class[i]>4:
        ax2.text(s=int(Class[i]), x=0.5,y=b+Class[i]/2-2.5, horizontalalignment='center',
                 verticalalignment='center', c='white',size=10,weight='bold')
    p=sorted_indivs[sorted_indivs[:,0]==i,2]
    x=sum(p==1)
    y=sum(p==2)
    z=sum(p==3)
    l=[x,y,z]
    c=0
    for j in range(3):
        if l[j]>4:
            ax2.text(s='%i'%l[j], x=2,y=b+c+l[j]/2-1, horizontalalignment='center',
            verticalalignment='center', c='white',size=10,weight='semibold')
        c+=l[j]
    b+=Class[i]
ax2.hlines(87.5,1.5,2.5,lw=0.5,colors='white')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

ax3=fig2.add_subplot(122)
sorted_indivs2=sorted(indivs,key=lambda x: (x[1], x[0]))

sorted_indivs2=np.array(sorted_indivs2)
for j in range(1,4):
    x=sorted_indivs2[sorted_indivs2[:,1]==j,0]
    print(CStat[j-1],'others',sum(x==0),'M1:',sum(x==1),'M2:',sum(x==2),'M3:',sum(x==3),'M4:',sum(x==4))
sorted_indivs2=np.vstack((sorted_indivs2.T,sorted_indivs2[:,1])).T
sorted_indivs2[:,1:3]+=4
g2=ax3.imshow(sorted_indivs2,cmap=cmaps,origin='lower',aspect='auto')
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel('Severity',size=13)
b=0
for i in range(0,3):
    ax3.text(s=CLabel[i], x=1.5,y=b+Severity[i]/2+0.9, horizontalalignment='center',
             verticalalignment='center', c='white',size=12,weight='bold')
    if i==0:
        ax3.text(s='%i'%(Severity[i]), x=1.5,y=b+Severity[i]/2-4.5, horizontalalignment='center',
                 verticalalignment='center', c='white',size=10,weight='bold')
    else:
        ax3.text(s='%i'%(Severity[i]), x=1.5,y=b+Severity[i]/2-2.5, horizontalalignment='center',
                 verticalalignment='center', c='white',size=10,weight='bold')
    p=sorted_indivs2[sorted_indivs2[:,1]==i+1,0]
    x=sum(p==0)
    y=sum(p==1)
    z=sum(p==2)
    l=[x,y,z]
    c=10
    for j in range(3):
        if l[j]>4:
            ax3.text(s='%i'%l[j], x=0,y=b+c+l[j]/2-1, horizontalalignment='center',
            verticalalignment='center', c='white',size=10,weight='semibold')
        c+=l[j]
    b+=Severity[i]
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
fig2.subplots_adjust(left=0.2,bottom=0.1,right=0.8,top=0.9,wspace=0.2,hspace=0.35)
fig2.savefig('Ratio1.svg')
fig2.savefig('Ratio1.png')
plt.show()