# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 09:25:39 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=160)

r=te.loada('''
model side_reaction
  J0: S + E -> SE; k1*k2*S*E - k2*ES;
  S = 5;
  E = 3;
  SE = E+S;
  k1 = 1.2;
  k2 = 0.4;
  ES = 1
end

model full_reaction
  var species S;
  A: side_reaction();
  B: side_reaction();
  A.S is S;
  B.S is S;
end
''')

r.draw()

b=te.loada('''
model side_reaction
  J0: S + E -> SE; k1*k2*S*E - k2*ES;
  S = 5;
  E = 3;
  SE = E+S;
  k1 = 1.2;
  k2 = 0.4;
  ES = 1
end''')
b.draw()
