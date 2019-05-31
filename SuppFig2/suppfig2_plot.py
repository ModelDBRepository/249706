
import numpy as np
import matplotlib.pylab as plt
import json

#######################################################

timeStep = 1        
deltaStim = 50
stimWindow = int(deltaStim/timeStep) 
Theta_eeMax_H = 25

datv = 'dat/'
dat_MD1 = open(datv+'supp_narrow.txt','r')
dat_MD2 = open(datv+'supp_broad.txt','r')
dt_MD1 = json.load(dat_MD1)
dt_MD2 = json.load(dat_MD2)
dat_MD1.close()
dat_MD2.close()


odi_4_narrow = np.array(dt_MD1[0])
odi_4_broad = np.array(dt_MD2[0])
odi_strt_narrow = np.array(dt_MD1[1])
odi_end_narrow = np.array(dt_MD1[2])
odi_strt_broad = np.array(dt_MD2[1])
odi_end_broad = np.array(dt_MD2[2])




bns=7
lw= 3
fs = 25


plt.figure()
plt.hist(odi_4_narrow,bins=bns,range=(-1,1),facecolor="None",edgecolor='b',label='Narrow',linewidth=lw)
plt.hist(odi_4_broad,bins=bns,range=(-1,1),facecolor="None",linestyle='dashed',edgecolor='r',label='Broad',linewidth=lw)
plt.legend(loc='best',fontsize=fs)
plt.xticks([-1,-0.5,0,0.5,1],fontsize=fs)
plt.yticks(fontsize=fs)
plt.ylim([0,1500])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.xlabel('ODI [a.u.]',fontsize=fs)
plt.ylabel('Number of neurons [a.u.]',fontsize=fs)
plt.title('Input ODI distribution',fontsize=fs)
plt.savefig('img/suppFig1_odi_distr'+'.eps',bbox_inches='tight')


plt.figure(figsize=(8,6))
plt.plot([0,1],[odi_strt_narrow[0],odi_end_narrow[0]],'-',marker='.',Markersize=25,linewidth=lw,color='b', label = 'Narrow')
plt.plot([0,1],[odi_strt_broad[0],odi_end_broad[0]],'-',marker='.',Markersize=25,linewidth=lw,color='r', label = 'Broad')
plt.xlim([-.3,1.3])
#plt.ylim([0,.65])
plt.xticks([0,1],['Before \n deprivation','After \n deprivation'],fontsize=fs)#,rotation=60
plt.ylabel('ODI',fontsize=fs)
plt.legend(loc='best',fontsize=int(0.7*fs))#'best')
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.title('Ocular dominance plasticity',fontsize=fs)
plt.savefig('img/suppFig1_odi_shift'+'.eps',bbox_inches='tight')



