import os
import numpy as np
import matplotlib.pylab as plt
#import matplotlib.lines as lns
import json
  

        
def load_obj(name ):
    with open('dat/' + name + '.txt', 'r') as f:
        return json.load(f)
    
def load_params(name ):
    script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    with open(os.path.join(script_dir,'params/' + name + '.txt'), 'r') as f:
        return json.load(f)
    

#####################################
###### Load
#####################################

m_par = load_params('model_params')

timeStep = m_par['timeStep']    
stimLen = m_par['stimLen']
deltaStim = m_par['deltaStim']
NrE_L3 = m_par['NrE_L3']
NrE_L4 = m_par['NrE_L4']
NrI_L3 = m_par['NrI_L3']
Nr4to3 = m_par['Nr4to3']
tau_e = m_par['tau_e']
tau_i = m_par['tau_i']
beta_v = m_par['beta_v']
beta_inh = m_par['beta_inh']
theta_v = m_par['theta_v']
vis_amp = m_par['vis_amp']
depr_value = m_par['depr_value']
nr_Ori = m_par['nr_Ori']

    
## Plasticity
p_par = load_params('plast_params')

wmax_taro = p_par['wmax_taro']
wmin_taro =  p_par['wmin_taro']
Theta_eeMax_H =  p_par['Theta_eeMax_H']
Theta_eeMax_L =  p_par['Theta_eeMax_L']
l_rate_ee =  p_par['l_rate_ee']
wmin_ei =  p_par['wmin_ei']
deltath =  p_par['deltath']
theta_ei_H =  p_par['theta_ei_H']
lrate_ei =  p_par['lrate_ei']
wmax_ie =  p_par['wmax_ie']
wmin_ie =  p_par['wmin_ie']


datvar = 'dat_totalSim_evo'
dt_w1 = load_obj(datvar)

L3_connections =   np.array(dt_w1['L3_connections']) # np.array(dt_w1[7])
w_ee_4to3 = np.array(dt_w1['w_ee_4to3']) 
w_ie_4to3 = np.array(dt_w1['w_ie_4to3']) 
w_ee_L3 = np.array(dt_w1['w_ee_L3']) 
w_ie_L3 = np.array(dt_w1['w_ie_L3']) 
w_ei_L3 = np.array(dt_w1['w_ei_L3']) 

#

####################################################################
bns=7
lw=3
fs = 25
lbl = 'Presynaptic Neurons'
lbl2 = 'Postsynaptic Neurons'

plt.figure()
dyn_fig = plt.imshow(w_ei_L3, interpolation='nearest', 
                            origin='bottom', 
                            aspect='auto', # get rid of this to have equal aspect
                            vmin=0,
                            vmax=np.amax(w_ei_L3), 
                            cmap='jet')
plt.xticks([0, 5,10,15],fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(lbl,fontsize=fs)
plt.ylabel(lbl2,fontsize=fs)
plt.title('Recurrent L II/III I-to-E weights',fontsize=fs)
dynfig_cb = plt.colorbar(ticks=[0.02, 0.2])
dynfig_cb.ax.set_yticklabels([0.02, 0.2],fontsize=fs) 
plt.savefig('img/NW_wei33'+'.eps',bbox_inches='tight')


plt.figure()
dyn_fig = plt.imshow(w_ie_L3, interpolation='nearest', 
                            origin='bottom', 
                            aspect='auto', # get rid of this to have equal aspect
                            vmin=wmin_ie,
                            vmax=wmax_ie, 
                            cmap='jet')
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(lbl,fontsize=fs)
plt.ylabel(lbl2,fontsize=fs)
plt.title('Recurrent L II/III E-to-I weights',fontsize=fs)
dynfig_cb = plt.colorbar(ticks=[wmin_ie+.01*wmax_ie, wmax_ie-.01*wmax_ie])
dynfig_cb.ax.set_yticklabels(['min w', 'max w'],fontsize=fs) 
plt.savefig('img/NW_wie33'+'.eps',bbox_inches='tight')


plt.figure()
dyn_fig = plt.imshow(w_ee_L3, interpolation='nearest', 
                            origin='bottom', 
                            aspect='auto', # get rid of this to have equal aspect
                            vmin=wmin_taro,
                            vmax=wmax_taro, 
                            cmap='jet')
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(lbl,fontsize=fs)
plt.ylabel(lbl2,fontsize=fs)
plt.title('Recurrent L II/III E-E weights',fontsize=fs)
dynfig_cb = plt.colorbar(ticks=[wmin_taro+.01*wmax_taro, wmax_taro-.01*wmax_taro])
dynfig_cb.ax.set_yticklabels(['min w', 'max w'],fontsize=fs) 
below = 3
for pp in range(5):
    ccols = [1-pp/5,pp/5,.2]
    line = plt.plot([2+pp*20,(pp+1)*20-2], [-below,-below], lw=10, color=ccols)
    line2 = plt.plot([-below,-below],[2+pp*20,(pp+1)*20-2], lw=10, color=ccols)
    line[0].set_clip_on(False)
    line2[0].set_clip_on(False)
plt.xlim([0,99])
plt.ylim([0,99])
ax = plt.gca()
ax.tick_params(axis='x', pad=11)
ax.tick_params(axis='y', pad=11)
plt.savefig('img/NW_wee33'+'.eps',bbox_inches='tight')

ee = w_ee_4to3[:,:]
ee[L3_connections==0] = -1
ee = np.ma.masked_where(ee < 0, ee)
cmap = plt.cm.jet
cmap.set_bad(color='white')
plt.figure()
dyn_fig2 = plt.imshow(ee, interpolation='nearest', 
                            origin='bottom', 
                            aspect='auto', # get rid of this to have equal aspect
                            vmin=wmin_taro,
                            vmax=wmax_taro, 
                            cmap=cmap)
plt.xticks(range(0,1000,300),fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(lbl,fontsize=fs)
plt.ylabel(lbl2,fontsize=fs)
plt.title('Feedforward E-E weights',fontsize=fs)
dynfig_cb2 = plt.colorbar(ticks=[wmin_taro+.01*wmax_taro, wmax_taro-.01*wmax_taro])
dynfig_cb2.ax.set_yticklabels(['min w', 'max w'],fontsize=fs) 
plt.savefig('img/NW_wee43'+'.eps',bbox_inches='tight')


L3_connections_inh = np.zeros((np.size(w_ie_4to3,axis=0),np.size(w_ie_4to3,axis=1)))# np.append(np.ones((NrI_L3,Nr4to3)),np.zeros((NrI_L3,NrE_L4-Nr4to3)),axis=1)
for ii in range(NrI_L3):
    L3_connections_inh[ii,:] =  np.tile(np.append(np.append(np.ones((1,int(.5*Nr4to3))),np.zeros((1,NrE_L4-Nr4to3))
                            ,axis=1),np.ones((1,int(.5*Nr4to3))),axis=1),(1,nr_Ori))
ie = w_ie_4to3[:,:]
ie[L3_connections_inh==0] = -1
ie = np.ma.masked_where(ie < 0, ie)
cmap = plt.cm.jet
cmap.set_bad(color='white')
plt.figure()
dyn_fig4 = plt.imshow(ie, interpolation='nearest', 
                            origin='bottom', 
                            aspect='auto', # get rid of this to have equal aspect
                            vmin=np.min(wmin_ie),
                            vmax=np.max(wmax_ie), 
                            cmap=cmap)
plt.xticks(range(0,1000,300),fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(lbl,fontsize=fs)
plt.ylabel(lbl2,fontsize=fs)
plt.xlim([0,1000])
plt.title('Feedforward E-to-I weights',fontsize=fs)
dynfig_cb4 = plt.colorbar(ticks=[wmin_ie+.01*wmax_ie, wmax_ie-.01*wmax_ie])
dynfig_cb4.ax.set_yticklabels(['min w', 'max w'],fontsize=fs) 
plt.savefig('img/NW_wie'+'.eps',bbox_inches='tight')


