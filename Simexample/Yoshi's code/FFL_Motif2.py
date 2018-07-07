# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 08:50:34 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=160)

#%% Repression FFL

print('\n')
print('*')*70
print('\n')

r=te.loada("""
    model FFL2_rep()
     M1 = 0;      M2 = 0;      M3 = 0;
     p_input = 0;      P2 = 0;      p_output = 0;
     L1 = .01;    L2 = .01;    L3 = .01;
     k1 = .65;    k2 = .65;    k3 = .65;
    dm1 = .5;    dm2 = .5;    dm3 = .5;
    dp_input = .5;    dp2 = .5;    dp_output = .5;
    TM1 = 15;    TM2 = 15;    TM3 = 15;
    Tr1 = .5;    Tr2 = .5;    Tr3 = .5;
     H1 =  1;     H2 =  1;     H3 =  2;

    R0: => M1   ; L1 + TM1 - dm1*M1
    R2: => M2   ; L2 + TM2 * 1/(1+(p_input/k2)^H2) - dm2*M2
    R3: => M3   ; L3 + TM3 * 1/(1+ (p_input/k3)^H3)* (P2/k3)^H3/(1+(P2/k3)^H3) - dm3*M3
    R4: => p_input   ; Tr1*M1 - dp_input*p_input
    R5: => P2   ; Tr2*M2 - dp2*P2
    R6: => p_output   ; Tr3*M3 - dp_output*p_output
  
    end
""")
r.reset()
#plt.close("all")
#res = r.simulate(0,50,1000)
#r.plot()
r.draw()
#r.reset()
r.exportToAntimony('FFL2_ant.txt') #export as antimony
