import sys
import numpy as np
#import matplotlib.pylab as plt
import json
from argparse import ArgumentParser
import os.path

#####################################
parser = ArgumentParser()
parser.add_argument("--mode",dest="mode",
                    help="Choose deprivation mode: MD-CL (=contra monocular deprivation), BD  (=binocular deprivation), MI  (=contra monocular inactivation), MD-IL  (=ipsi monocular deprivation).",
                    default='BD')
args = parser.parse_args()

#####################################
###### Functions
#####################################
    
def f_IEvol_WilCow(y,ipt,exc,iiw,ipw,iew,tauy):
    y = y + (timeStep/tauy)*(-y + gainF_relu( np.dot(ipw,ipt)+np.dot(iew,exc)-np.dot(iiw,y),beta_i ) )
    return y
    
def f_EEvol_WilCow(y,ipt,inh,eew,ipw,eiw,tauy):
    y = y + (timeStep/tauy)*(-y + gainF_relu( np.dot(ipw,ipt)+np.dot(eew,y)
                            -np.dot(eiw,inh),beta_e ) )
    return y
    
def gainF_relu(x,slopef):
    out = slopef*(x-theta_v)
    out = out*(x>theta_v)    
    return out
     
def f_excitatory_plasticity(pre=0,post=0,w=0,recurrent='no'):
    
    dw = ((np.outer(post,pre)>theta_ee_L)*(np.outer(post,pre)-theta_ee_H)
            *(np.outer(post,pre)-theta_ee_L))
            
    dw = np.tanh(dw/0.0001)        
        
    w = w + timeStep*l_rate_ee*dw  
    
    w[w>wmax_exc] = wmax_exc        
    w[w<wmin_exc] = wmin_exc
    
    if recurrent == 'yes':
        if np.size(w,axis=0) == np.size(w,axis=1):
            w = w - np.diag(np.diag(w))
        else:
            print('Problem with recurrent connections: no square matrix!')

    return w

def f_calc_odi(w_ee_L3b,w_ee_L4to3b,w_ei_L3b,w_ie_ffwb,w_ie_L3_effb,w_ii_L3b):
    # Store firing rates in order to calculate ocular dominance index
    tttime=50
    respv = np.zeros((NrE_L3,2))    # storage of data
    respv2 = np.zeros((NrE_L4,2))    # storage of data
    storV = np.zeros((NrE_L3,int(tttime/timeStep)))    # storage of data
    storVI = np.zeros((NrI_L3,int(tttime/timeStep)))    # storage of data
    storV4 = np.zeros((NrE_L4,int(tttime/timeStep)))    # storage of data
    
    for ippp in range(2):
        EvL4 = np.zeros(NrE_L4)          # vector of excitatory activities
        EvL3 = np.zeros(NrE_L3)          # vector of excitatory activities
        IvL3 = np.zeros(NrI_L3)          # vector of inhibitory activities
        
        if ippp == 1:
            # 'ipsi_closed'
            a = 0
            b = 1
        else:
            #'contra_closed'
            a = 1
            b = 0
             
        for tt in np.arange(0,tttime,timeStep):
            if tt%stimLen == 0:
                Iptv =  0*np.ones(NrE_L4)
                if tt%deltaStim == 0:
                    Visual_Ipt_gauss_ipsi = L4_connections[:,0]*np.array([a*vis_amp]*NrE_L4) 
                    Visual_Ipt_gauss_contra = L4_connections[:,1]*np.array([b*vis_amp]*NrE_L4)             
                    Visual_Ipt_gauss = a*Visual_Ipt_gauss_ipsi + b*Visual_Ipt_gauss_contra
                    Iptv = Visual_Ipt_gauss + Iptv 
                Iptv = 1*Iptv*(Iptv>0)
                  
            EvL4 = beta_e*(Iptv) 
            IvL3 = f_IEvol_WilCow(IvL3,layerfactor*EvL4,EvL3,w_ii_L3b,w_ie_ffwb,w_ie_L3_effb,tau_i)
            EvL3 = f_EEvol_WilCow(EvL3,layerfactor*EvL4,IvL3,w_ee_L3b,w_ee_L4to3b,w_ei_L3b,tau_e)  
            
            storV[:,int(tt/timeStep)] = EvL3
            storVI[:,int(tt/timeStep)] = IvL3
            storV4[:,int(tt/timeStep)] = EvL4
           
        ts = int(stimLen/timeStep)
        for xx in range(np.size(storV,axis=0)):
            respv[xx,ippp]= np.mean(storV[xx,:ts],axis=0)
        for xx in range(np.size(storV4,axis=0)):
            respv2[xx,ippp]= np.mean(storV4[xx,:ts],axis=0)
        
    return respv,respv2

#####################################
###### Parameters
#####################################

#*** General
timeStep = 1              # timestep [ms]
tmax = 10000                # Final time of simulation [ms]
#if args.mode == "MDi":
#    tmax = 5*tmax
storeDT = timeStep          # Frequency for storing data

stimLen = 25                # Duration of a visual input [ms]
deltaStim = 50              # Time between onset of two visual inputs [ms]

deprtime = 500              # Time of deprivation onset [ms]
RedInhibTime = 2000         # Time when inhibition will be reduced [ms]

#*** NEURONS
NrE_L3 = 1                  # total population of neurons in l3 
NrE_L4 = 1250                # total population of neurons in l4
NrExcNrns = 125              # nr of l4 -> l3 connections (per eye)
Nr4to3 = 2*NrExcNrns        # nr of l4 -> l3 connections (both eyes)
NrI_L3 = 1                  # total nr of inhibitory neurons in l3

tau_e = 1                   # excitatory time constant [ms]
tau_i = 1                   # inhibitory time constant [ms]
beta_e = .3                 # excitatory gain function slope
beta_i = .3                 # inhibitory gain function slope 
theta_v = 0                 # activity threshold 

#*** Activity
vis_amp = 10                # Amplitude of visual input
background_current = 10 #7.5    # Amplitude of background input

if args.mode == "MD-CL":
    Contra_depr =  0
    Ipsi_depr = 1
    bckg_depr = 1
elif args.mode == "BD":
    Contra_depr =  0
    Ipsi_depr = 0
    bckg_depr = 1
elif args.mode == "MI":
    Contra_depr =  0
    Ipsi_depr = 1
    bckg_depr = 0
elif args.mode == "MD-IL":
    Contra_depr = 1
    Ipsi_depr = 0
    bckg_depr = 1
else:
    sys.exit("ERROR: wrong mode argument, please choose between MD-CL, BD, MI or MD-IL")
    
theta_ee_H = 25 #13
theta_ee_L = 0 
layerfactor = 1             # Extra ffw input factor from layer 4 to 3


#*** Plasticity
l_rate_ee = 4e-6

#*** initial w
wmax_exc = 7.5/Nr4to3
wmin_exc = .05*wmax_exc

wee_ini = 0
wee_L4to3_ini = wmax_exc

wie_ini = 0

w_ie_ffw_ini = 7.5/Nr4to3  

wei_ini = 1.75
w_ii_ini = 0*1.5


#####################################
###### Variables
#####################################

#*** LAYER 3
Ev_L3 = np.zeros(NrE_L3)          # vector of excitatory activities
Iv_L3 = np.zeros(NrI_L3)          # vector of inhibitory activities

w_ie_ffw = w_ie_ffw_ini*np.ones((NrI_L3,NrE_L4))     # ffw weight matrix      
w_ei_L3 = wei_ini*np.ones((NrE_L3,NrI_L3))     # ffw weight matrix
w_ie_L3 = wie_ini*np.ones((NrI_L3,NrE_L3)) 
w_ii_L3 = w_ii_ini*np.ones((NrI_L3,NrI_L3))

w_ee_4to3 = wee_L4to3_ini*np.ones((NrE_L3,NrE_L4))
w_ee_L3 = wee_ini*np.ones((NrE_L3,NrE_L3))


#*** LAYER 4
Ev_L4 = np.zeros(NrE_L4)          # vector of excitatory activities
L4_connections = np.ones((NrE_L4,2))
L3_connections = np.zeros((NrE_L3,NrE_L4))

meancontra = 30
stdvcontra = 25

q = (meancontra +np.random.randn(NrE_L4)*stdvcontra)
q[q<0] = 100*np.random.rand(len(q[q<0]))
q[q>100] = 100*np.random.rand(len(q[q>100])) #100
q = np.sort(q)
q = np.expand_dims(q,1)


if os.path.exists('dat/connections.txt'):    
    dat_load = open('dat/connections.txt','r')
    dat_conncs = json.load(dat_load)
    dat_load.close()
    L3_connections = np.array(dat_conncs[0])
    L4_connections = np.array(dat_conncs[1])
else:
    L4_connections = np.append(q,100-q,axis=1)*.01
    arr = np.array((np.random.permutation(np.arange(1250))<250)*1)    
    L3_connections[0,:] = arr
    
L3_connections_inh =  np.append(np.append(np.ones((1,int(.6*Nr4to3))),np.zeros((1,NrE_L4-Nr4to3))
                            ,axis=1),np.ones((1,int(.4*Nr4to3))),axis=1)



storeV = np.zeros((15,int(tmax/storeDT)))    # storage of data
storeW = np.zeros((Nr4to3,int(tmax/storeDT)))    # storage of data


######################################
#----- Run -----#
######################################

stimWindow = int(deltaStim/timeStep) 
strV3 = np.zeros((NrE_L3,stimWindow))    # storage of data
strV4 = np.zeros((NrE_L4,stimWindow))    # storage of data
strI3 = np.zeros((NrI_L3,stimWindow))

bckgr = background_current
deprived = 1
BD_depr = 1
for tt in np.arange(0,int(tmax),timeStep):
    if tt%stimLen == 0:
                Iptv =  np.ones(NrE_L4)*0
                if tt%deltaStim == 0:
                    Visual_Ipt_gauss_ipsi = L4_connections[:,0]*np.array([BD_depr*vis_amp]*NrE_L4) 
                    Visual_Ipt_gauss_contra = L4_connections[:,1]*np.array([deprived*vis_amp]*NrE_L4)             
                    Visual_Ipt_gauss = Visual_Ipt_gauss_ipsi + Visual_Ipt_gauss_contra
                    Iptv = Iptv+Visual_Ipt_gauss+bckgr
                Iptv = 1*Iptv*(Iptv>0)
    
    storeT = int((tt%deltaStim)/timeStep)              
    strV3[:,storeT] = Ev_L3
    strV4[:,storeT] = Ev_L4
    strI3[:,storeT] = Iv_L3
         
    ie_ff = L3_connections_inh*w_ie_ffw #+ L3_connections_inh_low*w_ie_ffw_low
    ee_ff = L3_connections*w_ee_4to3 #+ L3_connections_low*w_ee_4to3_low
    
    if np.abs(tt-deprtime) < timeStep:
        rvv_d0,rvv2 = f_calc_odi(w_ee_L3,ee_ff,w_ei_L3,ie_ff,w_ie_L3,w_ii_L3)
        deprived = Contra_depr
        BD_depr = Ipsi_depr #depr_value  
        
        bckgr = bckgr*bckg_depr
        
    elif np.abs(tt-RedInhibTime) < timeStep:
        
        if args.mode=="MD-IL":
            w_ie_ffw = w_ie_ffw*.75
        else:
            w_ie_ffw = w_ie_ffw*.33
        
        
    Ev_L4 = beta_e*(Iptv)         
    Iv_L3 = f_IEvol_WilCow(Iv_L3,layerfactor*Ev_L4,Ev_L3,w_ii_L3,ie_ff,w_ie_L3,tau_i)
    Ev_L3 = f_EEvol_WilCow(Ev_L3,layerfactor*Ev_L4,Iv_L3,w_ee_L3,ee_ff,w_ei_L3,tau_e)  
    
    w_ee_4to3 = f_excitatory_plasticity(pre=Ev_L4,post=Ev_L3,w=w_ee_4to3,recurrent='no')
       
    if tt%storeDT < timeStep:
        storeV[0,int(tt/storeDT)] = Ev_L3[0]
        storeV[1,int(tt/storeDT)] = w_ee_4to3[0,0]
        storeV[2,int(tt/storeDT)] = w_ee_4to3[0,-1]
        storeV[3,int(tt/storeDT)] = Iv_L3[0]
        storeV[4,int(tt/storeDT)] = w_ee_4to3[0,10]
        storeV[5,int(tt/storeDT)] = w_ee_4to3[0,20]
        storeV[6,int(tt/storeDT)] = w_ee_4to3[0,40]
        storeV[7,int(tt/storeDT)] = w_ee_4to3[0,int(NrE_L4*.5)]
        storeV[8,int(tt/storeDT)] = np.mean(Ev_L4[:])
        storeV[9,int(tt/storeDT)] = Ev_L3[0]*Ev_L4[0]
        storeV[10,int(tt/storeDT)] = Ev_L3[0]*Ev_L4[-1]
        
        
        storeW[:,int(tt/storeDT)] = w_ee_4to3[0,L3_connections[0,:]==1]
        
      ##################################
      # Display the progress
      ###################################
    if np.abs(tt-tmax*0.05*np.round(20*tt/tmax)) < 0.5*timeStep:
        print('> '+str(int(np.round(100*tt/tmax)))+'%') 
      
rvv_d3,rvv2 = f_calc_odi(w_ee_L3,ee_ff,w_ei_L3,ie_ff,w_ie_L3,w_ii_L3)

##################################
# Store data
###################################

w_evo = storeW[:,0:-1:500]
odi_d0=(rvv_d0[:,1]-rvv_d0[:,0])/(rvv_d0[:,0]+rvv_d0[:,1]) 
odi_d3=(rvv_d3[:,1]-rvv_d3[:,0])/(rvv_d3[:,0]+rvv_d3[:,1])  



dataw  = [w_evo.tolist(),storeV[9,:].tolist(),storeV[10,:].tolist(),odi_d0.tolist(),
          odi_d3.tolist()]
data3 = open("dat/"+args.mode+'.txt','w')
json.dump(dataw,data3)
data3.close()


if not os.path.exists('dat/connections.txt'):    
    data2  = [L3_connections.tolist(),L4_connections.tolist()]
    datfile = open('dat/connections.txt','w')
    json.dump(data2,datfile)
    datfile.close()

##################################
# Figures
###################################
#
#lw=3
#fs = 18
#plt.figure()
#plt.pcolor(w_evo, vmin=wmin_exc, vmax=wmax_exc)
#cbar = plt.colorbar(ticks = [wmin_exc,wmax_exc])
#cbar.ax.set_yticklabels(['min','max'],fontsize=fs)
#plt.xlabel('Time [a.u.]',fontsize=fs)
#plt.ylabel('Synapses [a.u.]',fontsize=fs)
#plt.title('Weight evolution, '+args.mode, fontsize=fs)
#plt.text(2, 125, r'$\star$', fontsize=1.5*fs)
#plt.text(10, 125, r'$\ddag$', fontsize=fs)
#plt.xticks(fontsize=fs)
#plt.yticks(fontsize=fs)
#ax = plt.gca()
#ax.yaxis.set_label_coords(-.15,.5)
#plt.annotate("",
#            xy=(-15, 0), xycoords='data',
#            xytext=(-15, 250), textcoords='data',
#            arrowprops=dict(arrowstyle="<->",
#                            connectionstyle="arc3",linewidth=lw),annotation_clip=False)
#plt.text(-20, -10, 'Contra', fontsize=fs)
#plt.text(-17, 255, 'Ipsi', fontsize=fs)
##plt.savefig('img/Fig1_W_evo_'+args.mode+'.eps',bbox_inches='tight')
#
#plt.figure()
#plt.plot(np.arange(0,int(tmax),timeStep),storeV[3,:],linewidth=4,label='i3')
#plt.plot(np.arange(0,int(tmax),timeStep),storeV[0,:],linewidth=4,label='e3')
#plt.plot(np.arange(0,int(tmax),timeStep),storeV[8,:],linewidth=4,label='e4')
#plt.legend()
#
#
###
##
#
#plt.figure()
#ax = plt.plot(np.arange(0,int(tmax),timeStep),storeV[1,:],'-',linewidth=3,label='closed-eye dominated',color='b')
#plt.plot(np.arange(0,int(tmax),timeStep),storeV[2,:],'--',linewidth=3,label='open-eye dominated',color='r')
#plt.plot(np.arange(0,int(tmax),timeStep),storeV[7,:],'--',linewidth=3,label='open-eye dominated2',color='r')
#plt.legend(loc='best')
##plt.ylim([.2*wmax_exc,wmax_exc+.1*wmax_exc])
#plt.yticks([.1*wmax_exc,wmax_exc],['.1*wmax','wmax'],fontsize=fs)
#plt.xlabel('Time [a.u.]',fontsize=fs)
#plt.title('Weight evolution', fontsize=fs)
#
#
##stimWindow = deltaStim
#range1 = np.arange(0,stimWindow)
#range2 = np.arange(12*stimWindow,13*stimWindow)
#range3 = np.arange(150*stimWindow,151*stimWindow)
##range4 = np.arange(190*stimWindow,199*stimWindow)
#
#CL = np.append(np.append(storeV[9,range1],storeV[9,range2],axis=0),storeV[9,range3],axis=0)#,storeV[9,range4],axis=0)#
#IL = np.append(np.append(storeV[10,range1],storeV[10,range2],axis=0),storeV[10,range3],axis=0)#,storeV[10,range4],axis=0)#
#
#a = beta_e
#b = beta_i
#w_ef = wmax_exc*Nr4to3
#w_if = w_ie_ffw_ini*Nr4to3
#w_ee = wee_ini
#w_ei = wei_ini
#w_ie = wie_ini
#F = a*(vis_amp+0*background_current + background_current)
#
#denom = (1 - a*w_ee + a*b*w_ei*w_ie)
#Theory_prepost = layerfactor*F**2*(a*w_ef - a*b*w_ei*w_if)/denom + 0*F*(a-a*b*w_ei)*background_current/denom
#
#plt.figure()
#plt.plot(np.arange(0,int(3*deltaStim),timeStep),CL,'-',linewidth=3,label='closed-eye dominated',color='b')
#plt.plot(np.arange(0,int(3*deltaStim),timeStep),IL,'--',linewidth=3,label='open-eye dominated',color='r')
#plt.plot([0,int(3*deltaStim)],[theta_ee_H,theta_ee_H],'-',linewidth=2,label=r'$\theta_H$',color='k')
#plt.plot([0,int(3*deltaStim)],[theta_ee_L,theta_ee_L],'--',linewidth=2,label=r'$\theta_L$',color='k')
#plt.plot([0,int(deltaStim)],[Theory_prepost,Theory_prepost],'--',linewidth=2,label=r'theory',color='m')
#plt.legend(loc='best')
#plt.ylabel(r'$\rho_{pre} \rho_{post}$ [Hz$^2$]',fontsize= fs)
#plt.xticks([.05*deltaStim,1.05*deltaStim,2.5*deltaStim],['Normal','MD','MD+reduced inhib.'],fontsize=fs)
##plt.xlabel('Time [a.u.]',fontsize=fs)
#plt.title('Pre * post firing rates', fontsize=fs)
##plt.savefig('img/PrePost_'+args.mode+'.eps',bbox_inches='tight')
#
#bns=7



#tit_var = 'Ocular dominance plasticity'
#    
#plt.figure()
#plt.plot([0,1],[odi_d0[0],odi_d3[0]],'-',linewidth=lw,color='k')
#plt.plot([0,1],[odi_d0[0],odi_d3[0]],'.',Markersize=25,linewidth=lw,color='k')
#plt.xlim([-.1,1.1])
##plt.ylim([-.05,.3])
#plt.xticks([0,1],['before','after'],rotation=60,fontsize=fs)
#plt.ylabel('ODI',fontsize=fs)
#plt.title(tit_var,fontsize=fs)
##plt.savefig('img/ODI_'+args.mode+'.eps',bbox_inches='tight')
#
#odi_tot_4 = (np.sum(rvv2[:,1])-np.sum(rvv2[:,0]))/(np.sum(rvv2[:,1])+np.sum(rvv2[:,0]))
#odi_4 = (rvv2[:,1]-rvv2[:,0])/(rvv2[:,1]+rvv2[:,0])
#bns=7
#lw=3
#
#plt.figure()
#plt.subplot(1,2,1)
#plt.hist(odi_4,bins=bns,range=(-1,1),facecolor="None",edgecolor='b',label='Before',linewidth=lw)
##plt.legend()
#plt.title('4')
#
#plt.subplot(1,2,2)
#plt.bar([1],[np.sum(rvv2[:,1])],label='contra',color='b')
#plt.bar([2],[np.sum(rvv2[:,0])],label='ipsi',color='r')
#plt.legend()
#plt.title('4')



