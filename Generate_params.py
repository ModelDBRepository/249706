
import json
  
#####################################
###### functions
#####################################
        
def save_obj(obj, name ):
    with open('params/'+ name + '.txt', 'w') as f:
        json.dump(obj, f)


#####################################
###### General model parameters
#####################################
        
store_model = 'model_params'

PopExcL3 = 100      # total population of neurons in l3
PopExcL4 = 200      # total population of neurons in l4
NrInhNrns = 20*int(PopExcL3/100)  #total nr of l3 inhibitory neurons
Nr4to3 = 50         # nr of l4 -> l3 connections

nr_Ori = 5

#####

model_params = {
        'timeStep' : 1,        # timestep
        'stimLen' : 10,         # visual input duration
        'deltaStim' : 50 ,       # time between consecutive visual inputs
        'NrE_L4' : PopExcL4,
        'NrE_L3' : PopExcL3,
        'NrI_L3' : NrInhNrns,  
        'Nr4to3' : Nr4to3, 
        'tau_e' : 1,            # excitatory neuron time constant
        'tau_i' : 1,            # inhibitory neuron time constant
        'beta_v' : 0.3,          # slope constant (excitatory)
        'beta_inh' : 0.3,        # slope constant (inhibitory)
        'theta_v' : 0,          # activation threshold
        'vis_amp' : 10,         # visual input amplitude
        'bck_current' : 10,         # visual input amplitude
        'depr_value' : 0,     # deprivation factor
        'nr_Ori' : nr_Ori,           # number of orientations
        'layerfactor' : 10
        }

save_obj(model_params, store_model )


#####################################
###### Plasticity parameters
#####################################
        
store_plast = 'plast_params'

lrate_factor = 1#0.4
# E to E
wmax_taro = 1/Nr4to3
wmin_taro = 0.0005*wmax_taro
l_rate_ee = .4*lrate_factor*1.5*.025*.04/Nr4to3

# I to E
wmin_ei = .06#*50/Nr4to3
lrate_ei = 0.5*lrate_factor*.3*.5*.1/NrInhNrns

# E to I
wmax_ie = 1/Nr4to3  
wmin_ie = .2*wmax_ie/2.5
lrate_ie_Max= lrate_factor*3*.025*10*(.000005*nr_Ori)/Nr4to3

## initial weights
wee_ini = .1*wmax_taro
wee_L4to3_ini = .35*wmax_taro
wie_ini =  wmin_ie
w_ie_4to3_ini = .5*wmax_ie
wei_ini = 1*wmin_ei

####

plast_params = {
        'wmax_taro' : wmax_taro,        # 
        'wmin_taro' : wmin_taro,        # 
        'l_rate_ee' : l_rate_ee,        # 
        'Theta_eeMax_H' : 26,        # 
        'Theta_eeMax_L' : 15,        # 
        'wmin_ei' : wmin_ei,        # 
        'deltath' : .75,#35,        # 
        'theta_ei_H' : 6,        # 
        'lrate_ei' : lrate_ei,        # 
        'wmax_ie' : wmax_ie,        # 
        'wmin_ie' : wmin_ie,        # 
        'lrate_ie_Max' : lrate_ie_Max,        # 
        'ie_target' : 6,        # 
        'I3_avg_ini' : 3,        # 
        'tau_avg' : 12500,        # 
        'wee_ini' :wee_ini ,        # 
        'wee_L4to3_ini' : wee_L4to3_ini,        # 
        'wie_ini' : wie_ini,        # 
        'w_ie_4to3_ini' :w_ie_4to3_ini ,        # 
        'wei_ini' : wei_ini,        # 
        }

save_obj(plast_params, store_plast )

print('Parameters stored!')
