import sys
import numpy as np
#import matplotlib.pylab as plt
import json
import os
from argparse import ArgumentParser

#####################################

parser = ArgumentParser()
parser.add_argument("--mode",dest="mode",
                    help="Choose deprivation mode: MD-CL (=contra monocular deprivation), BD  (=binocular deprivation), MI  (=contra monocular inactivation), MD-IL  (=ipsi monocular deprivation).",
                    default='MD-CL')
args = parser.parse_args()

#####################################

#####################################
  
def save_obj(obj, name ):
    with open('dat/'+ name , 'w') as f:
        json.dump(obj, f)
        
def load_obj(name ):
    with open('dat/' + name + '.txt', 'r') as f:
        return json.load(f)
        
def load_params(name ):
    script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    with open(os.path.join(script_dir,'params/' + name + '.txt'), 'r') as f:
        return json.load(f)
  
np.random.seed()  

print('===================================')
print('START '+ args.mode +' sim')
print('===================================')

#####################################
###### Parameters
#####################################
store_name = 'dat_totalSim_preCP_'

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
vis_amp = 5#m_par['vis_amp']
depr_value = m_par['depr_value']
nr_Ori = m_par['nr_Ori']
NrE_L4_tot =nr_Ori*NrE_L4
bck_current = 15#m_par['bck_current']
layerfactor = m_par['layerfactor']

    
## Plasticity
p_par = load_params('plast_params')

wmax_taro = p_par['wmax_taro']
wmin_taro =  p_par['wmin_taro']
Theta_eeMax_H =  p_par['Theta_eeMax_H']
Theta_eeMax_L =  p_par['Theta_eeMax_L']
l_rate_ee =  p_par['l_rate_ee']
wmin_ei =  p_par['wmin_ei'] # 
deltath =  p_par['deltath']
theta_ei_H =  p_par['theta_ei_H']
lrate_ei =  p_par['lrate_ei']
wmax_ie =  p_par['wmax_ie']
wmin_ie =  0.4*wmax_ie #p_par['wmin_ie']
lrate_ie_Max =  p_par['lrate_ie_Max']
ie_target =  p_par['ie_target']
tau_avg =  p_par['tau_avg']

tmax1 = 0
tmax2 = 50001 #+20000
tmax = tmax1 + tmax2

deprtime = tmax1 +1000

calc_odi_time = 4000
odi_index = 0

storeDT = timeStep
stimWindow = int(deltaStim/timeStep) 
nr_per_Ori = int(NrE_L3/nr_Ori)



if args.mode == "MD-CL":
    Contra_depr =  depr_value#+ .0000001
    Ipsi_depr = 1
    bckg_depr = 1
elif args.mode == "BD":
    Contra_depr =  depr_value
    Ipsi_depr = depr_value
    bckg_depr = 1
elif args.mode == "MI":
    Contra_depr =  depr_value
    Ipsi_depr = 1
    bckg_depr = 0
elif args.mode == "MD-IL":
    Contra_depr = 1
    Ipsi_depr = depr_value
    bckg_depr = 1
else:
    sys.exit("ERROR: wrong mode argument, please choose between MD-CL, BD, MI or MD-IL")
    
CL_mod = 1
IL_mod = 1

#####################################
###### Load
#####################################
namestr= 'dat_totalSim_Evo_preCP'
data = load_obj(namestr)

w_ee_4to3 = np.array(data['w_ee_4to3']) # np.array(dt_w1[4])  
w_ee_L3 =  np.array(data['w_ee_L3']) #  np.array(dt_w1[0])  
w_ie_4to3 =  np.array(data['w_ie_4to3']) #  np.array(dt_w1[5])    
w_ei_L3 = np.array(data['w_ei_L3']) #  np.array(dt_w1[1])  
w_ie_L3 =  np.array(data['w_ie_L3']) # np.array(dt_w1[2])  
w_ii_L3 =  np.array(data['w_ii_L3']) # np.array(dt_w1[3]) #.1*np.ones((40,40))#
activity_shift =  np.array(data['activity_shift']) # np.array(dt_w1[6]) #  np.ones(PopExcL3) #.35*np.random.randn(PopExcL3) # 
L3_connections =   np.array(data['L3_connections']) # np.array(dt_w1[7])
L3_connections_inh =   np.array(data['L3_connections_inh']) # np.array(dt_w1[7])
I3_avg =  np.array(data['I3_avg']) # np.array(dt_w1[10])
L4_connections =  np.array(data['L4_connections']) # np.array(dt_w1[11])

corr_act = np.outer(activity_shift,np.array(np.ones(NrE_L4_tot)))
Ev_L3 = np.zeros(NrE_L3)          # vector of excitatory activities
Iv_L3 = np.zeros(NrI_L3)          # vector of inhibitory activities
Ev_L4 = np.zeros(NrE_L4_tot) 

#################################################################################

storeV = np.zeros((15,int(tmax/storeDT)))    # storage of data
weightV = np.zeros((15,int(tmax/storeDT)))    # storage of data
weightV2 = np.zeros((NrI_L3,int(tmax/storeDT)))    # storage of data
weightV3 = np.zeros((15,int(tmax/storeDT)))    # storage of data


#####################################
###### functions
#####################################

def f_IEvol_WilCow2(y,ipt,exc,iiw,ipw,ip_corr,iew,tauy,ori,bckg):
    tot_inputs = np.dot(ipw*ip_corr,ipt) 
    y = y + (timeStep/tauy)*(-y + gainF_relu( tot_inputs+np.dot(iew,exc)-np.dot(iiw,y)+bckg,beta_inh ) )
    return y
    
def f_EEvol_WilCow2(y,ipt,inh,eew,ipw,ip_corr,eiw,tauy,ori,bckg):
    tot_inputs = np.dot(ipw*ip_corr,ipt) 
    y = y + (timeStep/tauy)*(-y + gainF_relu( tot_inputs+np.dot(eew,y)
                            -np.dot(eiw,inh)+bckg,beta_v ) )
    return y
    
def gainF_relu(x,slopef):
    out = slopef*(x-theta_v)
    out = out*(x>theta_v)    
    return out
     
def f_e_to_i_plasticity(x,y,yavg,weight): 
    
    theta_bcm = yavg**2/ie_target
    
    aux = np.outer(y-theta_bcm,2*(x>0))*(np.outer(y-theta_bcm,2)<0)+np.outer(y-theta_bcm,x*(x>3.2))*(np.outer(y-theta_bcm,2)>0)
    dw = np.multiply(np.outer(y,np.ones(np.shape(x))),aux)
    
    weight = weight + timeStep*lrate_ie_Max*dw 
    
    weight[weight>wmax_ie] = wmax_ie
    weight[weight<wmin_ie] = wmin_ie       
  
    return weight

def f_i_to_e_plasticity(x,y,weight,activityshift): 
    pre_o = np.ones(np.shape(x))    
    dw = np.outer(y>theta_ei_H*activityshift-deltath,pre_o)*np.outer(y-theta_ei_H*activityshift,x) 
         
    weight = weight + timeStep*lrate_ei*dw 
    weight[weight<wmin_ei] = wmin_ei   
    return weight

def f_ee_plasticity(pre=0,post=0,whebb=0,recurrent='no'):
    
    dw = ((np.outer(post,pre)>Theta_eeMax_L)*(np.outer(post,pre)-Theta_eeMax_H)
            *(np.outer(post,pre)-Theta_eeMax_L))
            
    dw = np.tanh(dw/0.0001)        
        
    whebb = whebb + timeStep*l_rate_ee*dw  
    
    whebb[whebb>wmax_taro] = wmax_taro        
    whebb[whebb<wmin_taro] = wmin_taro
    
    if recurrent == 'yes':
        if np.size(whebb,axis=0) == np.size(whebb,axis=1):
            whebb = whebb - np.diag(np.diag(whebb))
        else:
            print('Problem with recurrent connections: no square matrix!')
  
    return whebb   
    
#% Calculate OSI
def f_OSI(rates):
    nr_Ori = np.size(rates,axis=1)
    rt_var = np.sort(rates,axis=1)
    epsilon = 1e-10
    tot_rt = np.divide(1,np.sum(rates,axis=1)+epsilon)
    ori = np.arange(0,2*np.pi,2*np.pi/nr_Ori)
    R = np.dot((np.tile(np.expand_dims(tot_rt,axis=1),(1,nr_Ori))*rt_var),np.exp(1j*ori))
    return np.abs(R)
    
    
def f_calc_odi(w_ee_L3b,w_ee_L4to3b,w_ei_L3b,w_ie_L4to3b,w_ie_L3_effb,w_ii_L3b,bckg):
    tttime= deltaStim#nr_Ori*deltaStim
    respv = np.zeros((NrE_L3,2))    # storage of data
    respv2 = np.zeros((NrE_L4_tot,2))    # storage of data
    respv3 = np.zeros((NrI_L3,2))    # storage of data
    storV = np.zeros((NrE_L3,int(tttime/timeStep),nr_Ori))    # storage of data
    storVI = np.zeros((NrI_L3,int(tttime/timeStep),nr_Ori))    # storage of data
    storV4 = np.zeros((NrE_L4_tot,int(tttime/timeStep),nr_Ori))    # storage of data
    
    for ippp in range(2):
        if ippp == 1:
            eye = 'ipsi_closed'
        else:
            eye = 'contra_closed'
            
        if eye == 'ipsi_closed':
            a = 0
            b = 1
        else:
            a = 1
            b = 0
            
        for ori in np.arange(nr_Ori):
            EvL4 = np.zeros(NrE_L4_tot)          # vector of excitatory activities
            EvL3 = np.zeros(NrE_L3)          # vector of excitatory activities
            IvL3 = np.zeros(NrI_L3)          # vector of inhibitory activities
            for tt in np.arange(0,deltaStim,timeStep):
                if tt%stimLen == 0:
                    Iptv =  np.zeros(NrE_L4_tot)
                    if tt%deltaStim == 0:
                        Ori = ori
                        aux_IL = np.zeros((nr_Ori,NrE_L4))
                        aux_IL[Ori,:] = vis_amp#+bckgr
                        aux_IL = np.reshape(aux_IL,(NrE_L4_tot))
                        aux_CL = np.zeros((nr_Ori,NrE_L4))
                        aux_CL[Ori,:] = vis_amp#+bckgr
                        aux_CL = np.reshape(aux_CL,(NrE_L4_tot))
                        Visual_Ipt_gauss_ipsi = L4_connections[:,0]*aux_IL
                        Visual_Ipt_gauss_contra = L4_connections[:,1]*aux_CL          
                        Visual_Ipt_gauss = a*Visual_Ipt_gauss_ipsi + b*Visual_Ipt_gauss_contra
                        Iptv = Visual_Ipt_gauss + Iptv
                    Iptv = 1*Iptv*(Iptv>0)
                    
                EvL4 = beta_v*(Iptv)
                
                IvL3 = f_IEvol_WilCow2(IvL3,layerfactor*EvL4,EvL3,w_ii_L3b,w_ie_L4to3b,L3_connections_inh,
                               w_ie_L3_effb,tau_i,ori,0)
                                
                EvL3 = f_EEvol_WilCow2(EvL3,layerfactor*EvL4,IvL3,w_ee_L3b,
                                          w_ee_L4to3b,L3_connections*corr_act,w_ei_L3b,tau_e,ori,0) 
                
                storV[:,int(tt/timeStep),ori] = EvL3
                storVI[:,int(tt/timeStep),ori] = IvL3
                storV4[:,int(tt/timeStep),ori] = EvL4
        
        maxTime = np.amax(storV,axis=1)
        maxOri = np.argmax(maxTime,axis=1)
        maxTime4 = np.amax(storV4,axis=1)
        maxOri4 = np.argmax(maxTime4,axis=1)
        maxTimeI = np.amax(storVI,axis=1)
        meanOriI = np.mean(maxTimeI,axis=1)
        
        for xx in range(np.size(storV,axis=0)):
            respv[xx,ippp]= maxTime[xx,maxOri[xx]] 
        for xx in range(np.size(storV4,axis=0)):
            respv2[xx,ippp]= maxTime4[xx,maxOri4[xx]]
        for xx in range(np.size(storVI,axis=0)):
            respv3[xx,ippp]= meanOriI[xx] 
    return respv,respv2,respv3
    
def f_calc_osi(w_ee_L3b,w_ee_L4to3b,w_ei_L3b,w_ie_L4to3b,w_ie_L3_effb,w_ii_L3b,nrOri,bckgr):
    tttime= deltaStim
    storV = np.zeros((NrE_L3,int(tttime/timeStep)))    # storage of data
    storVI = np.zeros((NrI_L3,int(tttime/timeStep)))    # storage of data
    storV4 = np.zeros((NrE_L4_tot,int(tttime/timeStep)))    # storage of data
    
    store_Osi_E = np.zeros((NrE_L3,nrOri))
    store_Osi_I = np.zeros((NrI_L3,nrOri))
    
    for ori in np.arange(nrOri):
        EvL4 = np.zeros(NrE_L4_tot)          # vector of excitatory activities
        EvL3 = np.zeros(NrE_L3)          # vector of excitatory activities
        IvL3 = np.zeros(NrI_L3)          # vector of inhibitory activities
        for tt in np.arange(0,deltaStim,timeStep):
            if tt%stimLen == 0:
                Iptv =  np.zeros(NrE_L4_tot)
                if tt%deltaStim == 0:
                    Ori = ori
                    aux_IL = np.zeros((nr_Ori,NrE_L4))
                    aux_IL[Ori,:] = vis_amp#+bckgr
                    aux_IL = np.reshape(aux_IL,(NrE_L4_tot))
                    aux_CL = np.zeros((nr_Ori,NrE_L4))
                    aux_CL[Ori,:] = vis_amp#+bckgr
                    aux_CL = np.reshape(aux_CL,(NrE_L4_tot))
                    Visual_Ipt_gauss_ipsi = L4_connections[:,0]*aux_IL
                    Visual_Ipt_gauss_contra = L4_connections[:,1]*aux_CL          
                    Visual_Ipt_gauss = Visual_Ipt_gauss_ipsi + Visual_Ipt_gauss_contra
                    Iptv = Visual_Ipt_gauss + Iptv 
                Iptv = 1*Iptv*(Iptv>0)
                            
            EvL4 = beta_v*(Iptv) 
            
            IvL3 = f_IEvol_WilCow2(IvL3,layerfactor*EvL4,EvL3,w_ii_L3b,w_ie_L4to3b,L3_connections_inh,
                           w_ie_L3_effb,tau_i,ori,0)
            
            EvL3 = f_EEvol_WilCow2(EvL3,layerfactor*EvL4,IvL3,w_ee_L3b,
                                      w_ee_L4to3b,L3_connections*corr_act,w_ei_L3b,tau_e,ori,0)  
          
            storV[:,int(tt/timeStep)] = EvL3
            storVI[:,int(tt/timeStep)] = IvL3
            storV4[:,int(tt/timeStep)] = EvL4
            
        store_Osi_E[:,ori] = np.amax(storV,axis=1)
        store_Osi_I[:,ori] = np.amax(storVI,axis=1)
        
    return store_Osi_E,store_Osi_I     
    
storeW = np.zeros((Nr4to3,int(tmax/storeDT)))    # storage of data
storeW2 = np.zeros((Nr4to3,int(tmax/storeDT)))    # storage of data

mx1 = np.argmax(w_ee_4to3[0,:])
mx2 = np.argmax(w_ee_4to3[1,:])
mx1 = mx1//NrE_L4
mx2 = mx2//NrE_L4
               
######################################
#----- Run -----#
######################################
Ori = int(np.floor(np.random.rand(1)*nr_Ori)[0])
w_ie_L3_eff = w_ie_L3


strV3 = np.zeros((NrE_L3,stimWindow))    # storage of data
strV4 = np.zeros((NrE_L4_tot,stimWindow))    # storage of data
strI3 = np.zeros((NrI_L3,stimWindow))


OSI_evolution_E = np.zeros(20)    
OSI_evolution_I = np.zeros(20)   
osi_index = 0     

maxI3 = 0    
store_odi = {}         

bckg = bck_current
for tt in np.arange(0,int(tmax),timeStep):
    if tt%stimLen == 0:
            Iptv =  np.zeros(NrE_L4_tot)
            if tt%deltaStim == 0:
                Ori = int(np.floor(np.random.rand(1)*nr_Ori)[0])
                Ori_vect = np.zeros((nr_Ori,NrE_L4))
                Ori_vect[Ori,:] = 1
                Ori_vect = np.reshape(Ori_vect,(NrE_L4_tot))
                Visual_Ipt_gauss_ipsi = L4_connections[:,0]*Ori_vect*IL_mod*vis_amp 
                Visual_Ipt_gauss_contra = L4_connections[:,1]*Ori_vect*CL_mod*vis_amp      
                Visual_Ipt_gauss = Visual_Ipt_gauss_ipsi + Visual_Ipt_gauss_contra
                Iptv = Visual_Ipt_gauss + Iptv +bckg*Ori_vect
            Iptv = 1*Iptv*(Iptv>0)
            
    if np.abs(tt-deprtime) < timeStep:
        print("START "+args.mode)
        
        rvv_d0, rvv2,rvvI = f_calc_odi(w_ee_L3,w_ee_4to3,w_ei_L3,w_ie_4to3,w_ie_L3_eff,w_ii_L3,bck_current)
        CL_mod = Contra_depr
        IL_mod = Ipsi_depr 
        bckg = bck_current*bckg_depr
            
        ie_target = ie_target*.1
     
    elif  np.abs((tt-tmax1)%calc_odi_time) < timeStep:
        print("calc odi "+str(odi_index))
        store_odi[str(odi_index)] = f_calc_odi(w_ee_L3,w_ee_4to3,w_ei_L3,w_ie_4to3,w_ie_L3_eff,w_ii_L3,bck_current)
        odi_index += 1
            
    strV3[:,int((tt%deltaStim)/timeStep)] = Ev_L3
    strV4[:,int((tt%deltaStim)/timeStep)] = Ev_L4
    strI3[:,int((tt%deltaStim)/timeStep)] = Iv_L3
    
    I3_avg = I3_avg + (timeStep/tau_avg)*(Iv_L3-I3_avg)
    
    w_ie_L3_eff = w_ie_L3 #*np.outer(h_inh_L3,np.ones(np.shape(Ev_L3)))
#        
    
    #####################################################################################################   
    # LAYER 4
    #####################################################################################################    
     
    Ev_L4 = beta_v*(Iptv)           
    
    #####################################################################################################   
    # LAYER 3
    #####################################################################################################        
            
    Iv_L3 = f_IEvol_WilCow2(Iv_L3,layerfactor*Ev_L4,Ev_L3,w_ii_L3,w_ie_4to3,
                            L3_connections_inh,w_ie_L3_eff,tau_i,Ori,0)
    Ev_L3 = f_EEvol_WilCow2(Ev_L3,layerfactor*Ev_L4,Iv_L3,w_ee_L3,
                           w_ee_4to3,L3_connections*corr_act,w_ei_L3,tau_e,Ori,0) 
                           
    
    
    w_ie_4to3 = f_e_to_i_plasticity(Ev_L4,Iv_L3,I3_avg,w_ie_4to3)
    w_ie_L3 = f_e_to_i_plasticity(Ev_L3,Iv_L3,I3_avg,w_ie_L3) 
    
    w_ee_L3 = f_ee_plasticity(pre=Ev_L3,post=Ev_L3,whebb=w_ee_L3,recurrent='no')
    w_ee_4to3 = f_ee_plasticity(pre=Ev_L4,post=Ev_L3,whebb=w_ee_4to3,recurrent='no')
                    
    w_ei_L3 = f_i_to_e_plasticity(Iv_L3,Ev_L3,w_ei_L3,activity_shift)  
    
    #####################################################################################################   
    #####################################################################################################   
    
    if Iv_L3[0]>500:
        print('BREAK')
        break

    
    if tt%storeDT < timeStep:
        
        storeV[0,int(tt/storeDT)] = Ev_L4[0]*Ev_L3[0] #Ev_L3[0]#
        storeV[1,int(tt/storeDT)] = Ev_L4[0]*Ev_L3[20] #Ev_L3[1]
        storeV[2,int(tt/storeDT)] = Ev_L4[0]*Ev_L3[40] #Iv_L3[0]
        storeV[3,int(tt/storeDT)] = Ev_L4[0]*Ev_L3[60] #Ev_L3[15]
        storeV[4,int(tt/storeDT)] = Ev_L4[0]#Ev_L4[25]
        storeV[5,int(tt/storeDT)] = Ev_L3[20]
        storeV[6,int(tt/storeDT)] = Ev_L3[40]#Ev_L3[48]
        storeV[7,int(tt/storeDT)] = Ev_L3[60]#max3[0]*max4[1]*spec_Ori[0,1,Ori_prev]
        storeV[8,int(tt/storeDT)] = w_ei_L3[0,0]
        storeV[9,int(tt/storeDT)] = w_ee_L3[0,1]
        storeV[10,int(tt/storeDT)] = w_ee_L3[0,NrE_L3-1]
        storeV[11,int(tt/storeDT)] = Iv_L3[0] #Iv_L3[0]*Ev_L3[0]
        storeV[12,int(tt/storeDT)] = Iv_L3[1] #Ev_L3[0]*Ev_L3[1]
        storeV[13,int(tt/storeDT)] = I3_avg[0]**2/ie_target #I3_avg[1] #Iv_L3[19] #Ev_L3[0]*Ev_L4[1]
        storeV[14,int(tt/storeDT)] = Ev_L4[0] #Ev_L3[0]*Ev_L3[-1]
#        
        weightV[0,int(tt/storeDT)] = w_ee_4to3[0,1]
        weightV[1,int(tt/storeDT)] = w_ie_4to3[0,0]
        weightV[2,int(tt/storeDT)] = w_ii_L3[0,0]
        weightV[3,int(tt/storeDT)] = Iv_L3[4] #Iv_L3[0]*Ev_L4[0]
        
        weightV2[:,int(tt/storeDT)] = w_ei_L3[0,:]
        
        store1 = w_ee_4to3[0,L3_connections[0,:]==1]
        store2 = w_ee_4to3[1,L3_connections[1,:]==1]
        storeW[:,int(tt/storeDT)] = store1[mx1*Nr4to3:(mx1+1)*Nr4to3]
        storeW2[:,int(tt/storeDT)] = store2[mx2*Nr4to3:(mx2+1)*Nr4to3]
        
      ##################################
      # Display the progress
      ###################################
    if np.abs(tt-tmax*0.05*np.round(20*tt/tmax)) < 0.5*timeStep:
#        print(w_ee_L3[0,1])
        print('> '+str(int(np.round(100*tt/tmax)))+'%') 
    
    nrpoints = .05
    if np.abs(tt-timeStep-tmax*nrpoints*np.round((1/nrpoints)*tt/tmax)) < 0.5*timeStep:
        osi_E,osi_I = f_calc_osi(w_ee_L3,w_ee_4to3,w_ei_L3,w_ie_4to3,w_ie_L3_eff,w_ii_L3,nr_Ori,bck_current)
        osi_E = f_OSI(osi_E)
        osi_I = f_OSI(osi_I)
        OSI_evolution_E[osi_index] = np.mean(osi_E)
        OSI_evolution_I[osi_index] = np.mean(osi_I)
        osi_index = osi_index + 1


epsilon = 1e-8
odi_tot_e = []
odi_tot_i = []
nr_stored = 0
for x,y in store_odi.items():
    e3 = store_odi[x][0]
    i3 = store_odi[x][2]
    odi_tot_e += [np.sum(e3[:,1])]
    odi_tot_e += [np.sum(e3[:,0])]
    odi_tot_i += [np.sum(i3[:,1])]
    odi_tot_i += [np.sum(i3[:,0])]
    nr_stored += 1
    

rvv_d0 = store_odi['0'][0]
rvv_d1 = store_odi[str(int(nr_stored/3))][0]
rvv_d2 = store_odi[str(2*int(nr_stored/3))][0]
rvv_d3 = store_odi[str(int(nr_stored)-1)][0]
    
#
odi_d0=(rvv_d0[:,1]-rvv_d0[:,0])/(rvv_d0[:,0]+rvv_d0[:,1]+epsilon) 
odi_d1=(rvv_d1[:,1]-rvv_d1[:,0])/(rvv_d1[:,0]+rvv_d1[:,1]+epsilon)  
odi_d2=(rvv_d2[:,1]-rvv_d2[:,0])/(rvv_d2[:,0]+rvv_d2[:,1]+epsilon)  
odi_d3=(rvv_d3[:,1]-rvv_d3[:,0])/(rvv_d3[:,0]+rvv_d3[:,1]+epsilon)  
odi_totd0 = (np.sum(rvv_d0[:,1])-np.sum(rvv_d0[:,0]))/(np.sum(rvv_d0[:,1])+np.sum(rvv_d0[:,0])+epsilon)
odi_totd1 = (np.sum(rvv_d1[:,1])-np.sum(rvv_d1[:,0]))/(np.sum(rvv_d1[:,1])+np.sum(rvv_d1[:,0])+epsilon)
odi_totd2 = (np.sum(rvv_d2[:,1])-np.sum(rvv_d2[:,0]))/(np.sum(rvv_d2[:,1])+np.sum(rvv_d2[:,0])+epsilon)
odi_totd3 = (np.sum(rvv_d3[:,1])-np.sum(rvv_d3[:,0]))/(np.sum(rvv_d3[:,1])+np.sum(rvv_d3[:,0])+epsilon)


wvar = store_name+args.mode+'.json'
w_v = {
        'odi_tot_e' : odi_tot_e,        # 
        'odi_tot_i' : odi_tot_i,
        'w_evo' : storeW[:,0:-1:500].tolist(),
        'rvv_d0' : rvv_d0.tolist(),        # 
        'rvv_d1' : rvv_d1.tolist(),     # 
        'rvv_d2' : rvv_d2.tolist(),     # 
        'rvv_d3' : rvv_d3.tolist(),        # 
        'w_ei' : weightV2.tolist()
        }

save_obj(w_v, wvar)

