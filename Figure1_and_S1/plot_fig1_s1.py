
import numpy as np
import matplotlib.pylab as plt
import json

#######################################################
def f_plot_weight_evo(w_evo,name,wmax_exc,wmin_exc): 
    lw=3
    fs = 18
    plt.figure()
    plt.pcolor(w_evo, vmin=wmin_exc, vmax=wmax_exc,cmap='jet')
    cbar = plt.colorbar(ticks = [wmin_exc,wmax_exc])
    cbar.ax.set_yticklabels(['min w','max w'],fontsize=fs)
    plt.xlabel('Time [a.u.]',fontsize=fs)
    plt.ylabel('Synapses [a.u.]',fontsize=fs)
    plt.title('Weight evolution, '+name, fontsize=fs)
    plt.text(1, -21, r'$\star$', fontsize=2.2*fs)
    plt.text(3.5, -23, r'$\ddag$', fontsize=1.5*fs)
    plt.xticks([0,5,10,15,20],fontsize=fs)
    plt.yticks(fontsize=fs)
    ax = plt.gca()
    ax.yaxis.set_label_coords(-.15,.5)
    plt.annotate("",
                xy=(-5, 0), xycoords='data',
                xytext=(-5, 250), textcoords='data',
                arrowprops=dict(arrowstyle="<->",
                                connectionstyle="arc3",linewidth=lw),annotation_clip=False)
    plt.text(-8, -10, 'Contra', fontsize=fs)
    plt.text(-7, 255, 'Ipsi', fontsize=fs)
    plt.savefig('img/Fig1_W_evo_'+name+'.eps',bbox_inches='tight')
    

#######################################################

timeStep = 1        
deltaStim = 50
stimWindow = int(deltaStim/timeStep) 
Theta_eeMax_H = 25

datv = 'dat/'
dat_MI = open(datv+'MI.txt','r')
dat_BD = open(datv+'BD.txt','r')
dat_MD = open(datv+'MD-CL.txt','r')
dat_MDi = open(datv+'MD-IL.txt','r')
dt_MI = json.load(dat_MI)
dt_BD = json.load(dat_BD)
dt_MD = json.load(dat_MD)
dt_MDi = json.load(dat_MDi)
dat_MI.close()
dat_MD.close()
dat_MDi.close()
dat_BD.close()

w_evol = {'MD-IL': np.array(dt_MDi[0]),
          'MD-CL': np.array(dt_MD[0]),
          'BD': np.array(dt_BD[0]),
          'MI': np.array(dt_MI[0])
          }

prepost_contra_MD = np.array(dt_MD[1])
prepost_ipsi_MD = np.array(dt_MD[2])
prepost_contra_BD = np.array(dt_BD[1])
prepost_ipsi_BD = np.array(dt_BD[2])
prepost_contra_MI = np.array(dt_MI[1])
prepost_ipsi_MI = np.array(dt_MI[2])
prepost_contra_MDi = np.array(dt_MDi[1])
prepost_ipsi_MDi = np.array(dt_MDi[2])

odi_strt_MDi = np.array(dt_MDi[3])
odi_end_MDi = np.array(dt_MDi[4])
odi_strt_MD = np.array(dt_MD[3])
odi_end_MD = np.array(dt_MD[4])
odi_strt_BD = np.array(dt_BD[3])
odi_end_BD = np.array(dt_BD[4])
odi_strt_MI = np.array(dt_MI[3])
odi_end_MI = np.array(dt_MI[4])




##
lw=3
fs = 25

#*** initial w
wmax_exc = 7.5/250
wmin_exc = .05*wmax_exc

for k,v in w_evol.items():
    f_plot_weight_evo(v,k,wmax_exc,wmin_exc)


#
#time = np.arange(0,int(tmax),timeStep)
#stride = 2000
#rng = np.append(np.arange(0,int(tmax/timeStep),stride),[len(time)-1])
#rng2 = np.append(np.append([0],np.arange(int(stride/2),int(tmax/timeStep),stride)),[len(time)-1])


range1 = np.arange(stimWindow-1,2*stimWindow-1)
range2 = np.arange(11*stimWindow,12*stimWindow)
range3 = np.arange(50*stimWindow,51*stimWindow)
CL_MD = np.append(np.append(prepost_contra_MD[range1],prepost_contra_MD[range2],axis=0),prepost_contra_MD[range3],axis=0)
IL_MD = np.append(np.append(prepost_ipsi_MD[range1],prepost_ipsi_MD[range2],axis=0),prepost_ipsi_MD[range3],axis=0)
CL_BD = np.append(np.append(prepost_contra_BD[range1],prepost_contra_BD[range2],axis=0),prepost_contra_BD[range3],axis=0)
IL_BD = np.append(np.append(prepost_ipsi_BD[range1],prepost_ipsi_BD[range2],axis=0),prepost_ipsi_BD[range3],axis=0)
CL_MI = np.append(np.append(prepost_contra_MI[range1],prepost_contra_MI[range2],axis=0),prepost_contra_MI[range3],axis=0)
IL_MI = np.append(np.append(prepost_ipsi_MI[range1],prepost_ipsi_MI[range2],axis=0),prepost_ipsi_MI[range3],axis=0)
CL_MDi = np.append(np.append(prepost_contra_MDi[range1],prepost_contra_MDi[range2],axis=0),prepost_contra_MDi[range3],axis=0)
IL_MDi = np.append(np.append(prepost_ipsi_MDi[range1],prepost_ipsi_MDi[range2],axis=0),prepost_ipsi_MDi[range3],axis=0)

dlta = 5
CL_MD[stimWindow-dlta:stimWindow] = -1
CL_MD[2*stimWindow-dlta:2*stimWindow] = -1
CL_BD[stimWindow-dlta:stimWindow] = -1
CL_BD[2*stimWindow-dlta:2*stimWindow] = -1
CL_MI[stimWindow-dlta:stimWindow] = -1
CL_MI[2*stimWindow-dlta:2*stimWindow] = -1
CL_MDi[stimWindow-dlta:stimWindow] = -1
CL_MDi[2*stimWindow-dlta:2*stimWindow] = -1
IL_MD[stimWindow-dlta:stimWindow] = -1
IL_MD[2*stimWindow-dlta:2*stimWindow] = -1
IL_BD[stimWindow-dlta:stimWindow] = -1
IL_BD[2*stimWindow-dlta:2*stimWindow] = -1
IL_MI[stimWindow-dlta:stimWindow] = -1
IL_MI[2*stimWindow-dlta:2*stimWindow] = -1
IL_MDi[stimWindow-dlta:stimWindow] = -1
IL_MDi[2*stimWindow-dlta:2*stimWindow] = -1

ytcks = [0,15,30]
y_lim = [0,50]
title_pos = [35,40]
rot_lbl = 35

fig = plt.figure(figsize=(24,16))
fig.suptitle('Pre * post firing rates', fontsize=fs)
plt.subplot(2,2,1)
plt.plot(np.arange(0,int(3*deltaStim),timeStep),CL_MD,'-',linewidth=4,color='b')
plt.plot(np.arange(0,int(3*deltaStim),timeStep),IL_MD,'-',linewidth=3,color='r')
plt.plot([0,int(3*deltaStim)],[Theta_eeMax_H,Theta_eeMax_H],'--',linewidth=2,color='k')
#plt.plot([0,int(3*deltaStim)],[Theta_eeMax_L,Theta_eeMax_L],'--',linewidth=2,color='k')
plt.text(152.5, 25, r'$\theta_H$', fontsize=fs)
#plt.text(152.5, 12, r'$\theta_L$', fontsize=fs)
plt.ylim(y_lim)
plt.text(title_pos[0], title_pos[1], 'MD-CL', fontsize=fs)
#plt.text(152.5, 3.3, r'$\theta$', fontsize=fs)
plt.xticks([])
plt.yticks(ytcks,fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.subplot(2,2,2)
plt.plot(np.arange(0,int(3*deltaStim),timeStep),CL_BD,'-',linewidth=4,color='b')#,label='closed-eye dominated')
plt.plot(np.arange(0,int(3*deltaStim),timeStep),IL_BD,'-',linewidth=3,color='r')#,label='open-eye dominated')
plt.plot([0,int(1.5*deltaStim)],[Theta_eeMax_H,Theta_eeMax_H],'--',color='k')#,linewidth=2,label=r'$\theta$',color='k')
#plt.plot([0,int(3*deltaStim)],[Theta_eeMax_L,Theta_eeMax_L],'--',linewidth=2,color='k')
#plt.text(152.5, 25, r'$\theta_H$', fontsize=fs)
#plt.text(152.5, 12, r'$\theta_L$', fontsize=fs)
plt.plot([0,1],[-1,-1],'s',markersize=15,color='b',label='CL dominated')
plt.plot([0,1],[-1,-1],'s',markersize=15,color='r',label='IL dominated')
plt.xticks([])
plt.yticks([])
#plt.legend(bbox_to_anchor=(1.05, .75), loc=2, borderaxespad=0.,numpoints=1)
plt.legend(bbox_to_anchor=(1.1, 1.17),loc=1,borderaxespad=0., numpoints=1,fontsize=fs)
plt.ylim(y_lim)
plt.text(title_pos[0], title_pos[1], 'BD', fontsize=fs)
#plt.text(152.5, 3.3, r'$\theta$', fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.subplot(2,2,3)
plt.plot(np.arange(0,int(3*deltaStim),timeStep),CL_MI,'-',linewidth=4,color='b')
plt.plot(np.arange(0,int(3*deltaStim),timeStep),IL_MI,'-',linewidth=3,color='r')
plt.plot([0,int(3*deltaStim)],[Theta_eeMax_H,Theta_eeMax_H],'--',linewidth=2,color='k')
#plt.plot([0,int(3*deltaStim)],[Theta_eeMax_L,Theta_eeMax_L],'--',linewidth=2,color='k')
plt.text(152.5, 25, r'$\theta_H$', fontsize=fs)
#plt.text(152.5, 12, r'$\theta_L$', fontsize=fs)
plt.ylabel(r'$\rho_{pre} \rho_{post}$',fontsize= fs+5)
plt.text(-40, 80, r'[Hz$^2$]', fontsize=fs,rotation=90)
ax = plt.gca()
ax.yaxis.set_label_coords(-.13,1)
plt.ylim(y_lim)
plt.yticks(ytcks,fontsize=fs)
plt.text(title_pos[0], title_pos[1], 'MI', fontsize=fs)
#plt.text(152.5, 3.3, r'$\theta$', fontsize=fs)
plt.xticks([.2*deltaStim,1.2*deltaStim,2.5*deltaStim],['Normal',r'depriv. ($\star$)','depriv.\n+low inhib. ($\ddag$)'],rotation=rot_lbl,fontsize=fs)
ax = plt.gca()
ax.tick_params(axis='x', pad=7)
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.subplot(2,2,4)
plt.plot(np.arange(0,int(3*deltaStim),timeStep),CL_MDi,'-',linewidth=4,color='b')
plt.plot(np.arange(0,int(3*deltaStim),timeStep),IL_MDi,'-',linewidth=3,color='r')
plt.plot([0,int(3*deltaStim)],[Theta_eeMax_H,Theta_eeMax_H],'--',linewidth=2,color='k')
#plt.plot([0,int(3*deltaStim)],[Theta_eeMax_L,Theta_eeMax_L],'--',linewidth=2,color='k')
plt.text(152.5, 25, r'$\theta_H$', fontsize=fs)
#plt.text(152.5, 12, r'$\theta_L$', fontsize=fs)
#plt.legend(loc='best'
plt.ylim(y_lim)
plt.yticks([])
plt.text(title_pos[0], title_pos[1], 'MD-IL', fontsize=fs)
#plt.text(152.5, 3.3, r'$\theta$', fontsize=fs)
plt.xticks([.2*deltaStim,1.2*deltaStim,2.5*deltaStim],['Normal',r'depriv. ($\star$)','depriv.\n+low inhib. ($\ddag$)'],rotation=rot_lbl,fontsize=fs)
ax = plt.gca()
ax.tick_params(axis='x', pad=7)
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.savefig('img/Fig1_prepost'+'.eps',bbox_inches='tight')

bns=7
plt.figure(figsize=(8,6))
plt.plot([0,1],[odi_strt_MD[0],odi_end_MD[0]],'-',marker='.',Markersize=25,linewidth=lw,color='k', label = 'MD-CL')
plt.plot([0,1],[odi_strt_MDi[0],odi_end_MDi[0]],'-',marker='^',Markersize=15,linewidth=lw,color='k', label = 'MD-IL')
plt.plot([0,1],[odi_strt_MI[0],odi_end_MI[0]],':s',Markersize=20,markeredgecolor='k',markerfacecolor='w',markeredgewidth=lw,linewidth=lw,color='k', label = 'MI')
plt.plot([0,1],[odi_strt_BD[0],odi_end_BD[0]],'--s',Markersize=8,linewidth=lw,color='k' , label = 'BD')
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
plt.savefig('img/Fig1_odi'+'.eps',bbox_inches='tight')



