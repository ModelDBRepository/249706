
import numpy as np
import matplotlib.pylab as plt
import json
import os
  

        
def load_obj(name ):
    with open('dat/' + name , 'r') as f:
        return json.load(f)
        
def load_params(name ):
    script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    with open(os.path.join(script_dir,'params/' + name + '.txt'), 'r') as f:
        return json.load(f)
    

#####################################
###### Load
#####################################
depvar =   'MD-CL' 

for xtrastr in ['Adult','preCP']:
    
    extray = 35
        
    datvar3 = 'dat_totalSim_'+xtrastr+depvar+ '.json'
    
    dt_w1 = load_obj(datvar3)
    odi_tot_e = np.array(dt_w1['odi_tot_e'])
    odi_tot_i = np.array(dt_w1['odi_tot_i'])
    w_evo = np.array(dt_w1['w_evo'])
    rvv_d0 = np.array(dt_w1['rvv_d0'])
    rvv_d1 = np.array(dt_w1['rvv_d1'])
    rvv_d2 = np.array(dt_w1['rvv_d2'])
    rvv_d3 = np.array(dt_w1['rvv_d3'])
    
    p_par = load_params('plast_params')
    wmax_taro = p_par['wmax_taro']
    wmin_taro =  p_par['wmin_taro']
    
    ####################################################################
    bns=7
    lw=3
    fs = 25
    ####################################################################
    #
    #xloc = []
    #idx = 0
    #for i in range(1,int(len(odi_tot_e))+1):
    #    xloc += [idx+i%2]
    #    idx = idx + 1 + i%2
    #    
    #plt.figure(figsize=(16,8))
    #plt.suptitle("Evolution of population response",fontsize=fs)
    #plt.subplot(2,1,1)
    #bars1 = plt.bar(xloc,odi_tot_e,label='e',color=['r','b'])
    #plt.ylabel('Excitatory', fontsize=fs)
    #plt.yticks([0,100,200],fontsize=fs)
    #plt.xticks([])
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.subplot(2,1,2)
    #bars2 = plt.bar(xloc,odi_tot_i,label='i',color=['r','b'])
    #plt.ylabel('Inhibitory', fontsize=fs)
    #plt.yticks([0,100],fontsize=fs)
    #xtl = [str(x)+"%" for x in range(0,102,20)]
    #plt.xticks(np.arange(1.5,xloc[-1],int((xloc[-1]-2)/5)),xtl,fontsize=.7*fs)
    #plt.xlabel('Simulation time',fontsize=fs)
    #plt.text(9, 120, 'Contra', fontsize=fs,color='r')
    #plt.text(9, 100, 'Ipsi', fontsize=fs,color='b')
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.savefig('img/Fig6_'+xtrastr+'_ResponseEvo_'+depvar+'.eps',bbox_inches='tight')
    #
    #
    #####################################################################
    #
    #plt.figure()
    #plt.pcolor(w_evo, vmin=wmin_taro, vmax=wmax_taro, cmap='jet')
    #cbar = plt.colorbar(ticks = [wmin_taro,wmax_taro])
    #cbar.ax.set_yticklabels(['min','max'],fontsize=fs)
    #plt.xlabel('Time [a.u.]',fontsize=fs)
    #plt.ylabel('Synapses [a.u.]',fontsize=fs)
    #plt.title('Weight evolution, '+depvar, fontsize=fs)
    #plt.xticks(fontsize=fs)
    #plt.yticks(fontsize=fs)
    #ax = plt.gca()
    #ax.yaxis.set_label_coords(-.26,.5)
    #plt.annotate("",
    #            xy=(-19, 0), xycoords='data',
    #            xytext=(-19, 50), textcoords='data',
    #            arrowprops=dict(arrowstyle="<->",
    #                            connectionstyle="arc3",linewidth=lw),annotation_clip=False)
    #plt.text(-38, -4, 'Contra', fontsize=fs)
    #plt.text(-31, 50, 'Ipsi', fontsize=fs)
    #plt.savefig('img/Fig6_'+xtrastr+'_W_evo_'+depvar+'.eps',bbox_inches='tight')
    
    ####################################################################
     
    epsilon = 1e-7
    odi_d0=(rvv_d0[:,1]-rvv_d0[:,0])/(rvv_d0[:,0]+rvv_d0[:,1]+epsilon) 
    odi_d1=(rvv_d1[:,1]-rvv_d1[:,0])/(rvv_d1[:,0]+rvv_d1[:,1]+epsilon)  
    odi_d2=(rvv_d2[:,1]-rvv_d2[:,0])/(rvv_d2[:,0]+rvv_d2[:,1]+epsilon)  
    odi_d3=(rvv_d3[:,1]-rvv_d3[:,0])/(rvv_d3[:,0]+rvv_d3[:,1]+epsilon)  
    odi_totd0 = (np.sum(rvv_d0[:,1])-np.sum(rvv_d0[:,0]))/(np.sum(rvv_d0[:,1])+np.sum(rvv_d0[:,0])+epsilon)
    odi_totd1 = (np.sum(rvv_d1[:,1])-np.sum(rvv_d1[:,0]))/(np.sum(rvv_d1[:,1])+np.sum(rvv_d1[:,0])+epsilon)
    odi_totd2 = (np.sum(rvv_d2[:,1])-np.sum(rvv_d2[:,0]))/(np.sum(rvv_d2[:,1])+np.sum(rvv_d2[:,0])+epsilon)
    odi_totd3 = (np.sum(rvv_d3[:,1])-np.sum(rvv_d3[:,0]))/(np.sum(rvv_d3[:,1])+np.sum(rvv_d3[:,0])+epsilon)
    
    #
    #plt.figure()
    #plt.hist(odi_d0,bins=bns,range=(-1,1),facecolor="None",edgecolor='b',label='Before',linewidth=lw)
    #plt.hist(odi_d3,bins=bns,range=(-1,1),facecolor="None",linestyle='dashed',edgecolor='r',label='After',linewidth=lw)
    ##plt.legend(loc=6,fontsize=fs)
    #plt.xticks([-1,-0.5,0,0.5,1],fontsize=fs)
    #plt.yticks(fontsize=fs)
    #plt.ylim([0,80-.6*extray])
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.annotate("",
    #            xy=(-.8, 73-.6*extray), xycoords='data',
    #            xytext=(.8, 73-.6*extray), textcoords='data',
    #            arrowprops=dict(arrowstyle="<->",
    #                            connectionstyle="arc3",linewidth=lw))
    #plt.text(-.9, 65-.6*extray, 'Ipsi', fontsize=fs)
    #plt.text(.65, 65-.6*extray, 'Contra', fontsize=fs)
    #plt.text(-1, 53-.7*extray, 'Before', fontsize=fs,color='b')
    #plt.text(-1, 46-.7*extray, 'After', fontsize=fs,color='r')
    #plt.xlabel('ODI [a.u.]',fontsize=fs)
    #plt.ylabel('Number of neurons [a.u.]',fontsize=fs)
    #plt.title('ODI distribution, '+xtrastr+' '+depvar,fontsize=fs)
    #plt.savefig('img/Fig6_odi_distri_'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')
    #
    #
    #plt.figure()
    #for ii in range(len(odi_d0)):
    #    diffv = odi_d0[ii]-odi_d3[ii]
    #    colorv = 2*np.abs(diffv)
    #    colorv = 1-(colorv*(colorv<.95)+.95*(colorv>.95))
    ##    plt.plot([0,1,2,3],[odi_d0[ii],odi_d1[ii],odi_d2[ii],odi_d3[ii]],linewidth=1,color=[.25,1-colorv,.25])
    #    if np.abs(diffv)<.05:
    #        plt.plot([0,1],[odi_d0[ii],odi_d3[ii]],color=[.75,.75,.75],linewidth=lw)
    #    else:
    #        if diffv > 0:
    #            plt.plot([0,1],[odi_d0[ii],odi_d3[ii]],color=[.25*colorv,.15+.75*colorv,.25*colorv],linewidth=lw)
    #        else:
    #            plt.plot([0,1],[odi_d0[ii],odi_d3[ii]],color=[.15+.75*colorv,.25*colorv,.25*colorv],linewidth=lw)
    #plt.xticks([0,1],['Before \n '+depvar,'After \n '+depvar],rotation = 60,fontsize=fs)
    #plt.yticks(fontsize=fs)
    #plt.ylabel('ODI [a.u.]',fontsize=fs)
    #plt.xlim([-.1,1.1])
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.title('Individual ODI shifts, '+xtrastr+' '+depvar,fontsize=fs)
    #plt.savefig('img/Fig6_indi_odi'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')
    #
    #plt.figure()
    ##plt.plot([0,1,2,3],[odi_totd0,odi_totd1,odi_totd2,odi_totd3],linewidth=lw,color='k')
    ##plt.errorbar([0,1,2,3],[odi_totd0,odi_totd1,odi_totd2,odi_totd3],[np.std(odi_d0,axis=0),
    ##             np.std(odi_d1,axis=0),np.std(odi_d2,axis=0),np.std(odi_d3,axis=0)],linewidth=lw,color='k')
    ##plt.xticks([0,1,2,3],[ '0%','33%','66%','100%'],rotation = 0,fontsize=fs)
    #plt.plot([0,3],[odi_totd0,odi_totd3],linewidth=lw,color='k')
    #plt.errorbar([0,3],[odi_totd0,odi_totd3],[np.std(odi_d0,axis=0)/np.sqrt(np.size(odi_d0,axis=0)),
    #             np.std(odi_d3,axis=0)/np.sqrt(np.size(odi_d3,axis=0))],linewidth=lw,color='k')
    #plt.xticks([0,3],[ 'Before \n'+depvar,'After \n'+depvar],rotation = 60,fontsize=fs)
    #plt.yticks(fontsize=fs)
    #plt.ylabel('Population ODI [a.u.]',fontsize=fs)
    ##plt.xlabel('Simulation time',fontsize=fs)
    #plt.xlim([-.4,3.4])
    #plt.ylim([-.5,1])
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.title('Population ODI, '+xtrastr+' '+depvar,fontsize=fs)
    #plt.savefig('img/Fig6_pop_odi'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')
    #
    #
    #plt.figure()
    #plt.bar([1,4],[np.mean(rvv_d0[:,1]),np.mean(rvv_d3[:,1])],label='contra',color='k',linewidth=lw,align='edge')
    #plt.bar([2,5],[np.mean(rvv_d0[:,0]),np.mean(rvv_d3[:,0])],label='ipsi',color='w',linewidth=lw,align='edge',edgecolor='k')
    ##plt.bar([1,2],[np.sum(rvv_d0[:,1]),np.sum(rvv_d3[:,1])],label='contra',color='b')
    ##plt.bar([3,4],[np.sum(rvv_d0[:,0]),np.sum(rvv_d3[:,0])],label='ipsi',color='r')
    #plt.legend(loc=9,fontsize=fs)
    #plt.xticks([2,5],['Before '+depvar, 'After '+depvar],fontsize=fs)
    #plt.xlim([.9,6.])
    #plt.yticks(fontsize=fs)
    #plt.ylabel('Mean rate [Hz]',fontsize=fs)
    #ax = plt.gca()
    #ax.spines['right'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    #plt.title('Population response, '+xtrastr+' '+depvar,fontsize=fs)
    #plt.savefig('img/Fig6_pop_response'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')
    #
    
    #######################################################
    
    datvar3 = 'suppdat_totalSim_'+xtrastr+depvar+ '.json'
    
    dt_w1 = load_obj(datvar3)
    odi_tot_e = np.array(dt_w1['odi_tot_e'])
    odi_tot_i = np.array(dt_w1['odi_tot_i'])
    w_evo = np.array(dt_w1['w_evo'])
    rvv_d0 = np.array(dt_w1['rvv_d0'])
    rvv_d1 = np.array(dt_w1['rvv_d1'])
    rvv_d2 = np.array(dt_w1['rvv_d2'])
    rvv_d3 = np.array(dt_w1['rvv_d3'])
    
    epsilon = 1e-7
    
    
    
    plt.figure() 
    
    plt.errorbar([0,3],[odi_totd0,odi_totd3],[np.std(odi_d0,axis=0)/np.sqrt(np.size(odi_d0,axis=0)),
                 np.std(odi_d3,axis=0)/np.sqrt(np.size(odi_d3,axis=0))],linewidth=lw,color='k')
    
    odi_d0=(rvv_d0[:,1]-rvv_d0[:,0])/(rvv_d0[:,0]+rvv_d0[:,1]+epsilon) 
    odi_d1=(rvv_d1[:,1]-rvv_d1[:,0])/(rvv_d1[:,0]+rvv_d1[:,1]+epsilon)  
    odi_d2=(rvv_d2[:,1]-rvv_d2[:,0])/(rvv_d2[:,0]+rvv_d2[:,1]+epsilon)  
    odi_d3=(rvv_d3[:,1]-rvv_d3[:,0])/(rvv_d3[:,0]+rvv_d3[:,1]+epsilon)  
    odi_totd0 = (np.sum(rvv_d0[:,1])-np.sum(rvv_d0[:,0]))/(np.sum(rvv_d0[:,1])+np.sum(rvv_d0[:,0])+epsilon)
    odi_totd1 = (np.sum(rvv_d1[:,1])-np.sum(rvv_d1[:,0]))/(np.sum(rvv_d1[:,1])+np.sum(rvv_d1[:,0])+epsilon)
    odi_totd2 = (np.sum(rvv_d2[:,1])-np.sum(rvv_d2[:,0]))/(np.sum(rvv_d2[:,1])+np.sum(rvv_d2[:,0])+epsilon)
    odi_totd3 = (np.sum(rvv_d3[:,1])-np.sum(rvv_d3[:,0]))/(np.sum(rvv_d3[:,1])+np.sum(rvv_d3[:,0])+epsilon)
    
    
    
    #plt.plot([0,3],[odi_totd0,odi_totd3],linewidth=lw,color='k')
    plt.errorbar([0,3],[odi_totd0,odi_totd3],[np.std(odi_d0,axis=0)/np.sqrt(np.size(odi_d0,axis=0)),
                 np.std(odi_d3,axis=0)/np.sqrt(np.size(odi_d3,axis=0))],linestyle='--',linewidth=lw,color='k')
    plt.xticks([0,3],[ 'Before \n'+depvar,'After \n'+depvar],rotation = 0,fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.ylabel('Population ODI [a.u.]',fontsize=fs)
    #plt.xlabel('Simulation time',fontsize=fs)
    plt.xlim([-.4,3.4])
    plt.ylim([-.5,1])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    plt.title('Population ODI, '+xtrastr+' '+depvar,fontsize=fs)
    plt.savefig('img/Fig6_pop_odi_both'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')
    
    plt.figure()
    plt.hist(odi_d0,bins=bns,range=(-1,1),facecolor="None",edgecolor='b',label='Before',linewidth=lw)
    plt.hist(odi_d3,bins=bns,range=(-1,1),facecolor="None",linestyle='dashed',edgecolor='r',label='After',linewidth=lw)
    #plt.legend(loc=6,fontsize=fs)
    plt.xticks([-1,-0.5,0,0.5,1],fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.ylim([0,80-.6*extray])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    plt.annotate("",
                xy=(-.8, 73-.6*extray), xycoords='data',
                xytext=(.8, 73-.6*extray), textcoords='data',
                arrowprops=dict(arrowstyle="<->",
                                connectionstyle="arc3",linewidth=lw))
    plt.text(-.9, 65-.6*extray, 'Ipsi', fontsize=fs)
    plt.text(.65, 65-.6*extray, 'Contra', fontsize=fs)
    plt.text(-1, 53-.7*extray, 'Before', fontsize=fs,color='b')
    plt.text(-1, 46-.7*extray, 'After', fontsize=fs,color='r')
    plt.xlabel('ODI [a.u.]',fontsize=fs)
    plt.ylabel('Number of neurons [a.u.]',fontsize=fs)
    plt.title('ODI distribution, '+xtrastr+' '+depvar,fontsize=fs)
    plt.savefig('img/suppFig6_odi_distri_'+xtrastr+'_'+depvar+'.eps',bbox_inches='tight')



