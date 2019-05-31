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
                    default='MD-IL')
args = parser.parse_args()

#####################################
  
def save_obj(obj, name ):
    with open('dat/'+ name , 'w') as f:
        json.dump(obj, f)
        
def load_obj(name ):
    with open('dat/' + name, 'r') as f:
        return json.load(f)
        
def load_params(name ):
    script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    with open(os.path.join(script_dir,'params/' + name + '.txt'), 'r') as f:
        return json.load(f)
  
np.random.seed()  

nr_trials = 50
store_odi_diff = np.zeros(nr_trials)
d0 = [[] for _ in range(nr_trials)]
d3 = [[] for _ in range(nr_trials)]
acs = [[] for _ in range(nr_trials)]


for trialnr in range(nr_trials):
    
    print('===================================')
    print(str(trialnr+1)+' of '+str(nr_trials))
    print('===================================')
        
    #####################################
    ###### Parameters
    #####################################
    wvar = 'dat_totalSim_thresholds_'
    
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
    NrE_L4_tot =nr_Ori*NrE_L4
    bck_current = m_par['bck_current']
    layerfactor = m_par['layerfactor']
    
        
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
    lrate_ie_Max =  p_par['lrate_ie_Max']
    ie_target =  p_par['ie_target']
    tau_avg =  p_par['tau_avg']
    I3avg_init = p_par['I3_avg_ini']
    wee_ini = p_par['wee_ini'] 
    wee_L4to3_ini = p_par['wee_L4to3_ini'] 
    wie_ini = p_par['wie_ini'] 
    w_ie_4to3_ini = p_par['w_ie_4to3_ini'] 
    wei_ini = p_par['wei_ini'] 
          
    
    Theta_eeMax_H =  Theta_eeMax_H + np.random.randn(1)*0.05*Theta_eeMax_H
    Theta_eeMax_L =  Theta_eeMax_L + np.random.randn(1)*0.05*Theta_eeMax_L
    
    tmax1 = 50000
    tmax2 = 50001 
    tmax = tmax1 + tmax2
    
    dayvar = int((tmax2-1000)/300) #1 'day' of deprivation
    deprtime = tmax1 +1000
    
    calc_odi_time = 5000
    odi_index = 0
    
    storeDT = timeStep
    stimWindow = int(deltaStim/timeStep) 
    nr_per_Ori = int(NrE_L3/nr_Ori)
    
    
    if args.mode == "MD-CL":
        Contra_depr =  depr_value
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
    I3_avg = I3avg_init*np.ones(NrI_L3)
    
    activity_shift = 1 + .1*np.random.randn(NrE_L3) 
    corr_act = np.outer(activity_shift,np.array(np.ones(NrE_L4_tot)))
    
    backge = .1#1#0
    backgi = .1#1
    unos = (np.diag(np.ones(nr_Ori)) + backge*(np.ones((nr_Ori,nr_Ori))-np.diag(np.ones(nr_Ori)))) 
    unos[0,-1] = backge
    unos[-1,0] = backge
    unos=np.expand_dims(unos,axis=1)
    spec_Ori = np.repeat(np.repeat(unos,nr_per_Ori,axis=0),NrE_L4,axis=1)
    spec_Ori = spec_Ori + np.random.rand(nr_Ori*nr_per_Ori,NrE_L4,nr_Ori)*.5
    
    if NrI_L3>1:
        unos = backgi*(np.ones((NrI_L3,nr_Ori)))
        for ii in np.arange(NrI_L3):
            unos[ii,int(np.floor(np.random.rand()*nr_Ori))] = 1
    #    print(unos)
        unos=np.expand_dims(unos,axis=1)
        spec_Ori_inh = np.repeat(unos,NrE_L4,axis=1)
    else:
        spec_Ori_inh = np.ones((NrI_L3,NrE_L4,nr_Ori))
    
    
    
    #*** Other
    NrIp = 2       # total nr of eyes
        
    #####################################
    ###### Variables
    #####################################
    
    #################################################################################
    # LAYER 3
    #################################################################################
    Ev_L3 = np.zeros(NrE_L3)          # vector of excitatory activities
    Iv_L3 = np.zeros(NrI_L3)          # vector of inhibitory activities
    
    w_ie_4to3 = w_ie_4to3_ini*spec_Ori_inh#*np.ones((NrI_L3,NrE_L4,nr_Ori))     # #ffw weight matrix   
    w_ie_4to3[w_ie_4to3<wmin_ie] = wmin_ie 
    w_ie_4to3 =np.reshape(w_ie_4to3,(NrI_L3,nr_Ori*NrE_L4),order='F')
    w_ei_L3 = wei_ini*np.ones((NrE_L3,NrI_L3))     # ffw weight matrix
    w_ie_L3 = wie_ini*np.ones((NrI_L3,NrE_L3)) 
    w_ii_L3 = np.zeros((NrI_L3,NrI_L3))
    
    L3_connections = np.zeros((NrE_L3,NrE_L4_tot))
    
    subfam = int(NrE_L4/5)
    
    for ii in range(0,int(.1*NrE_L3)):
        arr = []
        distr = [5,5,5,10,25]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    for ii in range(int(.1*NrE_L3),int(.2*NrE_L3)):
        arr = []
        distr = [5,10,10,10,15]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    for ii in range(int(.2*NrE_L3),int(.4*NrE_L3)):
        arr = []
        distr = [10,10,10,10,10]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    for ii in range(int(.4*NrE_L3),int(.6*NrE_L3)):
        arr = []
        distr = [10,15,10,10,5]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    for ii in range(int(.6*NrE_L3),int(.8*NrE_L3)):
        arr = []
        distr = [10,15,15,5,5]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    for ii in range(int(.8*NrE_L3),int(1*NrE_L3)):
        arr = []
        distr = [25,10,5,5,5]
        distr = np.array(distr)*Nr4to3/50
        sp = len(distr)
        for _ in range(nr_Ori):
            for jj in range(sp):
                a = (np.random.permutation(np.arange(subfam))<distr[jj])*1
                arr = arr+a.tolist()
        L3_connections[ii,:] = np.array(arr)
    
    np.random.shuffle(L3_connections)
        
    L3_connections_inh = np.zeros((NrI_L3,NrE_L4_tot))# np.append(np.ones((NrI_L3,Nr4to3)),np.zeros((NrI_L3,NrE_L4-Nr4to3)),axis=1)
    for ii in range(NrI_L3):
        L3_connections_inh[ii,:] =  np.tile(np.append(np.append(np.ones((1,int(.6*Nr4to3))),np.zeros((1,NrE_L4-Nr4to3))
                                ,axis=1),np.ones((1,int(.4*Nr4to3))),axis=1),(1,nr_Ori))
    
    w_ee_4to3 = wee_L4to3_ini*spec_Ori#np.ones((NrE_L3,NrE_L4))
    w_ee_4to3[w_ee_4to3<wmin_taro] = wmin_taro#np.ones((NrE_L3,NrE_L4))
    w_ee_4to3 = np.reshape(w_ee_4to3,(NrE_L3,nr_Ori*NrE_L4),order='F')
    w_taroee_4to3 = w_ee_4to3
    w_ee_L3 = wee_ini+wee_ini*np.random.randn(NrE_L3,NrE_L3) #np.ones((NrE_L3,NrE_L3))
    w_ee_L3[w_ee_L3<wmin_taro] = wmin_taro
    
    #################################################################################
    
    
    #################################################################################
    # LAYER 4
    #################################################################################
    Ev_L4 = np.zeros(NrE_L4_tot)          # vector of excitatory activities
    L4_connections = []
    meancontra = 30
    stdvcontra = 35
    for _ in range(nr_Ori):
        q = meancontra +np.random.randn(NrE_L4)*stdvcontra
        q[q<0] = 100*np.random.rand(len(q[q<0])) #0
        q[q>100] = 100
        q = (np.sort(q)).tolist()
        L4_connections.append(q)
    q = np.reshape(L4_connections,[NrE_L4_tot])
    L4_connections = .01*np.transpose(np.array([q,100-q]))
    
    
    #################################################################################
    
    storeV = np.zeros((15,int(tmax/storeDT)))    # storage of data
    weightV = np.zeros((15,int(tmax/storeDT)))    # storage of data
    weightV2 = np.zeros((15,int(tmax/storeDT)))    # storage of data
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
            
            wmin_ie = .4*wmax_ie # see methods
            
            rvv_d0, rvv2,rvvI = f_calc_odi(w_ee_L3,w_ee_4to3,w_ei_L3,w_ie_4to3,w_ie_L3_eff,w_ii_L3,bck_current)
            CL_mod = Contra_depr
            IL_mod = Ipsi_depr 
            bckg = bck_current*bckg_depr
                
            ie_target = ie_target*.1
                
        elif  np.abs((tt-tmax1)%calc_odi_time) < timeStep and tt>tmax1:
            print("calc odi "+str(odi_index))
            store_odi[str(odi_index)] = f_calc_odi(w_ee_L3,w_ee_4to3,w_ei_L3,w_ie_4to3,w_ie_L3_eff,w_ii_L3,bck_current)
            odi_index += 1
                
        strV3[:,int((tt%deltaStim)/timeStep)] = Ev_L3
        strV4[:,int((tt%deltaStim)/timeStep)] = Ev_L4
        strI3[:,int((tt%deltaStim)/timeStep)] = Iv_L3
        
        I3_avg = I3_avg + (timeStep/tau_avg)*(Iv_L3-I3_avg)
        
        w_ie_L3_eff = w_ie_L3
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
    bns=7
    lw=3
    
    store_odi_diff[trialnr] = odi_totd3-odi_totd0
    d0[trialnr] = odi_d0.tolist()
    d3[trialnr] = odi_d3.tolist()
    acs[trialnr] = activity_shift.tolist()


    w_v = {
            'diff_ODI' : store_odi_diff.tolist()     ,
            'rvv_d0' : d0,        # 
            'rvv_d3' : d3,        # 
            'activity_shift': acs,
            'trialnr': trialnr
                }
    
    save_obj(w_v, wvar+args.mode+'_ODIdiff_'+str(nr_trials) +'.txt')

