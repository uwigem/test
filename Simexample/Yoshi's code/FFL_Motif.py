# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:00:49 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=160)

#%% Feed-Forward Loop Motif

print('\n')
print('*')*70
print('\n')

r = te.loada('''
//create a FFL: p_input->P2->p_output ; p_input->P4->p_output
//assume all activation TF 
//assume time units are in minutes

model ffl1()

    #reactions
    
    trnsc_m1: $g1 => m1 ; a_m1 * g1 + L1
    deg_m1: m1 => ; d_m1*m1
    input: m1 => m1 + p_input ; a_p_input*m1 //input located in tandem with protein production (entry)
    deg_p_input: p_input => ; d_p_input*p_input 
    
    trnsc_m2: $g2 + p_input => p_input + m2 ; v1 * g2 + L2
    deg_m2: m2 => ; d_m2*m2
    trnsl_p2: m2 => m2 + p2 ; a_p2*m2
    deg_p2: p2 => ; d_p2*p2
    
    trnsc_m3: $g3 + p2 + p_input => p2 + p_input + m3 ; v2 * g3 + L3
    deg_m3: m3 => ; d_m3*m3
    trnsl_p_output: m3 => m3 + p_output ; a_p_output*m3
    output: p_output => ; d_p_output*p_output //output located in tandem with protein degradation (exit)
    
    #reaction rates
    v1 := (p_input^h1)/(K1^h1 + p_input^h1)
    v2 := Vm2 * (p2^h2)/(K2^h2 + p2^h2) * (Vm3 * (p_input^h3)/(K3^h3 + p_input^h3))
    
    #species and parameters
    g1 = 1 ; g2 = 1  ; g3 = 1  ; //DNA count
    m1 = 0 ; m2 = 0  ; m3 = 0  ; //mRNA count
    p_input = 0 ; p2 = 0  ; p_output = 0  ; //protein count
    
    #a_m, d_m from bionumbers DNA binding model
    
    L1 = .1     ; L2 = .1     ; L3 = .1     ; //leak rate
    a_m1 = 20   ; a_m2 = 20   ; a_m3 = 20   ; //transc rate
    a_p_input = 1    ; a_p2 = 1    ; a_p_output = 1    ; //transl rate
    d_m1 = .6   ; d_m2 = .6   ; d_m3 = .6   ; //mRNA deg rate
    d_p_input = 0.09 ; d_p2 = 0.09 ; d_p_output = 0.09 ; //prot deg rate
    
    Vm1 = 20 ; Vm2 = 20 ; Vm3 = 20 ; //max transc rate
    h1 = 2   ; h2 = 2   ; h3 = 2   ; //hill coeff
    K1 = .2  ; K2 = .2  ; K3 = .2  ; //dissociation constant
end
''')


r.reset()
#plt.close("all")
#res = r.simulate(0,50,1000)
#r.plot()
#r.draw()
#r.reset()
r.exportToAntimony('FFL1_ant.txt') #export as antimony

#%% matplotlib plotting
#time = res[:,0]
#Pro2 = res[:,2]
#Pro3 = res[:,4]
#Pro4 = res[:,6]
#
#plt.figure(1)
#plt.plot(time,Pro2)
#plt.plot(time,Pro3)
#plt.plot(time,Pro4)
#plt.legend(['Pro2','Pro3','Pro4'])