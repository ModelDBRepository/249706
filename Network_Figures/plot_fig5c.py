
import numpy as np
import matplotlib.pylab as plt
import json
import os

        
def load_obj(name ):
    with open('dat/' + name, 'r') as f:
        return json.load(f)
        
def load_params(name ):
    script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    with open(os.path.join(script_dir,'params/' + name + '.txt'), 'r') as f:
        return json.load(f)
    

#####################################
###### Load
#####################################


datvar3 = 'dat_totalSim_MD-CL.json'
dt_w1 = load_obj(datvar3)
datvar3 = 'dat_totalSim_MD-IL.json'
od_IL = load_obj(datvar3)
datvar3 = 'dat_totalSim_MI.json'
od_MI = load_obj(datvar3)
datvar3 = 'dat_totalSim_BD.json'
od_BD = load_obj(datvar3)

rvv_d0 = np.array(dt_w1['rvv_d0'])
rvv_d3 = np.array(dt_w1['rvv_d3'])
d0_IL = np.array(od_IL['rvv_d0'])
d3_IL = np.array(od_IL['rvv_d3'])
d0_MI = np.array(od_MI['rvv_d0'])
d3_MI = np.array(od_MI['rvv_d3'])
d0_BD = np.array(od_BD['rvv_d0'])
d3_BD = np.array(od_BD['rvv_d3'])

epsilon = 1e-8
odi_totd0 = (np.sum(rvv_d0[:,1])-np.sum(rvv_d0[:,0]))/(np.sum(rvv_d0[:,1])+np.sum(rvv_d0[:,0]) + epsilon)
odi_totd3 = (np.sum(rvv_d3[:,1])-np.sum(rvv_d3[:,0]))/(np.sum(rvv_d3[:,1])+np.sum(rvv_d3[:,0]) + epsilon)
t0_IL = (np.sum(d0_IL[:,1])-np.sum(d0_IL[:,0]))/(np.sum(d0_IL[:,1])+np.sum(d0_IL[:,0]) + epsilon)
t3_IL = (np.sum(d3_IL[:,1])-np.sum(d3_IL[:,0]))/(np.sum(d3_IL[:,1])+np.sum(d3_IL[:,0]) + epsilon)
t0_MI = (np.sum(d0_MI[:,1])-np.sum(d0_MI[:,0]))/(np.sum(d0_MI[:,1])+np.sum(d0_MI[:,0]) + epsilon)
t3_MI = (np.sum(d3_MI[:,1])-np.sum(d3_MI[:,0]))/(np.sum(d3_MI[:,1])+np.sum(d3_MI[:,0]) + epsilon)
t0_BD = (np.sum(d0_BD[:,1])-np.sum(d0_BD[:,0]))/(np.sum(d0_BD[:,1])+np.sum(d0_BD[:,0]) + epsilon)
t3_BD = (np.sum(d3_BD[:,1])-np.sum(d3_BD[:,0]))/(np.sum(d3_BD[:,1])+np.sum(d3_BD[:,0]) + epsilon)

odi_d0=(rvv_d0[:,1]-rvv_d0[:,0])/(rvv_d0[:,0]+rvv_d0[:,1] + epsilon) 
odi_d3=(rvv_d3[:,1]-rvv_d3[:,0])/(rvv_d3[:,0]+rvv_d3[:,1] + epsilon)  
s0_IL=(d0_IL[:,1]-d0_IL[:,0])/(d0_IL[:,0]+d0_IL[:,1] + epsilon) 
s3_IL=(d3_IL[:,1]-d3_IL[:,0])/(d3_IL[:,0]+d3_IL[:,1] + epsilon)  
s0_MI=(d0_MI[:,1]-d0_MI[:,0])/(d0_MI[:,0]+d0_MI[:,1] + epsilon) 
s3_MI=(d3_MI[:,1]-d3_MI[:,0])/(d3_MI[:,0]+d3_MI[:,1] + epsilon)  
s0_BD=(d0_BD[:,1]-d0_BD[:,0])/(d0_BD[:,0]+d0_BD[:,1] + epsilon) 
s3_BD=(d3_BD[:,1]-d3_BD[:,0])/(d3_BD[:,0]+d3_BD[:,1] + epsilon)  


lw=3
fs = 25


plt.figure()
plt.plot([0,2],[odi_totd0,odi_totd3],'-',marker='.',Markersize=25,linewidth=lw,color='k', label = 'MD-CL')
plt.plot([0,2],[t0_IL,t3_IL],'-',marker='^',Markersize=15,linewidth=lw,color='k', label = 'MD-IL')
plt.plot([0,2],[t0_BD,t3_BD],':s',Markersize=20,markeredgecolor='k',markerfacecolor='w',
         markeredgewidth=lw,linewidth=lw,color='k', label = 'MI')
plt.plot([0,2],[t0_MI,t3_MI],'--s',Markersize=8,linewidth=lw,color='k' , label = 'BD')


#plt.errorbar([0,2],[odi_totd0,odi_totd3],[np.std(odi_d0,axis=0)/np.sqrt(np.size(odi_d0,axis=0)),
#             np.std(odi_d3,axis=0)/np.sqrt(np.size(odi_d3,axis=0))],linestyle='-',marker='.',Markersize=25,linewidth=lw,color='k', label = 'MD-CL')
#plt.errorbar([0,2],[t0_IL,t3_IL],[np.std(s0_IL,axis=0)/np.sqrt(np.size(s0_IL,axis=0)),
#             np.std(s3_IL,axis=0)/np.sqrt(np.size(s3_IL,axis=0))],linewidth=lw,color='r')
#plt.errorbar([0,2],[t0_MI,t3_MI],[np.std(s0_MI,axis=0)/np.sqrt(np.size(s0_MI,axis=0)),
#             np.std(s3_MI,axis=0)/np.sqrt(np.size(s3_MI,axis=0))],linewidth=lw,color='r')
#plt.errorbar([0,2],[t0_BD,t3_BD],[np.std(s0_BD,axis=0)/np.sqrt(np.size(s0_BD,axis=0)),
#             np.std(s3_BD,axis=0)/np.sqrt(np.size(s3_BD,axis=0))],linewidth=lw,color='r')
    
plt.xticks([0,2],[ 'Before','After'],rotation = 0,fontsize=fs)
plt.yticks(fontsize=fs)
plt.ylabel('Population ODI [a.u.]',fontsize=fs)
plt.legend(loc='best',fontsize=int(0.7*fs))
plt.xlim([-.4,3.4])
plt.ylim([-.5,1])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.title('Population ODI',fontsize=fs)
plt.savefig('img/Fig4_pop_odi_all.eps',bbox_inches='tight')


