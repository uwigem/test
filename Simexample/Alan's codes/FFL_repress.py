# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:16:49 2018

@author: Alan
"""

import tellurium as te
import roadrunner
import matplotlib.pyplot as plt
import numpy as np
#   Activation    a2 + (P1/k2)/(1+P1/k2)^H2
#   Repression    a3 + 1/(1+ (P1/k3)^H3)
r=te.loada("""
    model FFL_repress()
     M1 = 0;      M2 = 0;      M3 = 0;
     P1 = 0;      P2 = 0;      P3 = 0;
     L1 = .01;    L2 = .01;    L3 = .01;
     k1 = .65;    k2 = .65;    k3 = .65;
    dm1 = .5;    dm2 = .5;    dm3 = .5;
    dp1 = .5;    dp2 = .5;    dp3 = .5;
    TM1 = 15;    TM2 = 15;    TM3 = 15;
    Tr1 = .5;    Tr2 = .5;    Tr3 = .5;
     H1 =  1;     H2 =  1;     H3 =  2;

    R0: => M1   ; L1 + TM1 - dm1*M1
    R2: => M2   ; L2 + TM2 * 1/(1+(P1/k2)^H2) - dm2*M2
    R3: => M3   ; L3 + TM3 * 1/(1+ (P1/k3)^H3)* (P2/k3)^H3/(1+(P2/k3)^H3) - dm3*M3
    R4: => P1   ; Tr1*M1 - dp1*P1
    R5: => P2   ; Tr2*M2 - dp2*P2
    R6: => P3   ; Tr3*M3 - dp3*P3
  
    end
""")
tmax=200
result = r.simulate(0, tmax, 200,)
plt.figure()
plt.grid(color='k', linestyle='-', linewidth=.4)
plt.ylim(0,np.max(result[:,4:7])*1.1)
plt.xlim(0,tmax)
plt.yticks(np.arange(0,np.max(result[:,4:7])*1.1,np.max(result[:,4:7])/12))

#M1 , = plt.plot (result[:,0],result[:,1], label = 'M1')
#M2 , = plt.plot (result[:,0],result[:,3], label = 'M2')
#M2 , = plt.plot (result[:,0],result[:,6], label = 'M3')
P1 , = plt.plot (result[:,0],result[:,4], label = 'P1')
P2 , = plt.plot (result[:,0],result[:,5], label = 'P2')
P3 , = plt.plot (result[:,0],result[:,6], label = 'P3')
plt.legend([P1, P2, P3], ['P1', 'P2', 'P3'])

#plt.legend([M1, M2, M3, P1, P2, P3], ['M1', 'M2', 'M3', 'P1', 'P2', 'P3'])