
import numpy as np
import matplotlib.pylab as plt
import json
#import time
#from IPython import display
  
def save_obj(obj, name ):
    with open('dat/'+ name + '.txt', 'w') as f:
        json.dump(obj, f)
        
def load_obj(name):
    with open('dat/' + name + '.txt', 'r') as f:
        return json.load(f)
  
    
plot_mode = ['MD-CL','BD', 'MI', 'MD-IL']#,'MD-CL','MD-CL']#, 'BD', 'MI', 'MD-IL']
nrsims = 50

fs = 25
lw = 3

data = []
for name in plot_mode:
    a = load_obj('dat_totalSim_thresholds_'+name+'_ODIdiff_'+str(nrsims))
    newdat = [(-1*np.array(a['diff_ODI'])).tolist()]
#    print(newdat[0][0:15])
    data = data+[newdat[0][:39]]
    

plt.figure()
plt.boxplot(data,whis=1.5, boxprops={'linewidth':lw}, medianprops={'linewidth':lw}, 
            whiskerprops={'linewidth':lw}, capprops={'linewidth':lw}, flierprops={'markerfacecolor':'k','markeredgewidth':lw})
plt.xticks(range(1,1+len(plot_mode)),plot_mode,fontsize=fs)
plt.yticks(fontsize=fs)
plt.ylabel('ODI shift', fontsize=fs)
plt.title('Robustness for threshold choice', fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
plt.savefig('img/Fig_thresholds.eps',bbox_inches='tight')




a = load_obj('dat_totalSim_thresholds_MD-CL_ODIdiff_'+str(nrsims))
newdat_1 = np.reshape(np.array(a['rvv_d0']),[-1])
#print(len(newdat1))
newdat_2 = np.reshape(np.array(a['rvv_d3']),[-1])
newdat_3 = np.reshape(np.array(a['activity_shift']),[-1])


newdat1 = newdat_1#[kk*50:ll*50]
newdat2 = newdat_2#[kk*50:ll*50]
newdat3 = newdat_3#[kk*50:ll*50]

diffv = newdat2 - newdat1
act_ratios = newdat3

act_pos = act_ratios[diffv>.15]


plt.figure()
plt.boxplot([act_pos,act_ratios],whis=1.5, boxprops={'linewidth':lw}, medianprops={'linewidth':lw}, 
        whiskerprops={'linewidth':lw}, capprops={'linewidth':lw}, flierprops={'markerfacecolor':'k','markeredgewidth':lw})
plt.plot([0.6,1.9],[1,1],'--',linewidth=2,color='r')
plt.text(1.1,1.01,'Pop. mean',color='r',fontsize=fs)
plt.xticks([1,2],['ctr-int. \n shifters', 'population'],rotation = 45, fontsize = fs)
plt.yticks([.8,1,1.2],fontsize=fs)
plt.title('Positively shifting neurons \n fire at low rates',fontsize=fs)
plt.ylabel('Relative firing rate',fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
plt.savefig('img/Counter_shifts_pop'+'.eps',bbox_inches='tight')



plt.figure()
plt.plot(act_ratios,diffv,'.',color = 'k', markersize=5)
plt.xlabel('Ratio of neuron rate to \n mean population rate',fontsize=fs)
plt.ylabel(r'$\Delta$ ODI [a.u.]',fontsize=fs)
plt.title('Positively shifting neurons \n fire at low rates',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
##

diff_ipsi = diffv[newdat1<-.75]
diff_contra = diffv[newdat1>-.75]
diffdat = [diff_ipsi, diff_contra]

plt.figure()
plt.boxplot(diffdat,whis=1.5, boxprops={'linewidth':lw}, medianprops={'linewidth':lw}, 
        whiskerprops={'linewidth':lw}, capprops={'linewidth':lw}, flierprops={'markerfacecolor':'k','markeredgewidth':lw})
plt.xticks([1,2],['<-0.75','>-0.75'],fontsize=fs)
plt.yticks([-1,0,1],fontsize=fs)
plt.title('OD shift vs. original ODI',fontsize=fs)
plt.ylabel('ODI shift [a.u.]',fontsize=fs)
plt.xlabel('Original ODI [a.u.]',fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
plt.savefig('img/ODIvsOrig_pop'+'.eps',bbox_inches='tight')


plt.figure()
plt.plot(newdat1,diffv,'.',color = 'k', markersize=5)
plt.xlabel('Original ODI [a.u.]',fontsize=fs)
plt.ylabel(r'$\Delta$ ODI [a.u.]',fontsize=fs)
plt.title('OD shift vs. original ODI',fontsize=fs)
plt.xticks([-1,-.5,0,.5,1],fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')


