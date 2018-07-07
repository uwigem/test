# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np
import antimony
import string
import itertools

np.set_printoptions(linewidth=160)

def product_gen(n):
    for r in itertools.count(1):
        for i in itertools.product(n, repeat=r):
            yield "".join(i)

#%% Test model: putting two models together
numReads = 2
reads  = []
for i in np.arange(0,numReads):
    reads.append(te.readFromFile('FFL1_ant.txt'))
    reads.append(te.readFromFile('FFL2_ant.txt'))
    #combine this with for loop below this later
# Note: The string object name has to match the model name

readsalpha = zip(string.ascii_uppercase[:numReads],reads)

models = ""
species = ""
count  = 0
for i in np.arange(0,numReads):
    models = models + reads[i][1]
    species = species + "var species p_c"+str(count)+";\n"
    
    currStr = reads[i][1]
    modules = modules + str(reads[i][0])+ " : " + currStr.split("model",1)[1]
    
    count += count
    

combined =  models + '''
model combined
    var species p_c;
    #!!module names have to be identical to model names from imported models!!
    A : ffl1();
    B : ffl1();
    A.p_input is p_c;
    B.p_output is p_c;
end
'''

#What I want
#    combined = add(all models I have) + '''
#    model combined
#       specify a global input node
#    
#        A : model1();
#        B : model2();
#        repeat for all models I have
#        var species p_c;
#        repeat for n-1 modules
#        A.p_input is p_c;
#        B.p_output is p_c;
#        repeat for n-1 modules
#    end 
#    '''


#properly flatten the combined model
antimony.loadAntimonyString(combined)
sbmlstr = antimony.getSBMLString('combined')
antimony.loadSBMLString(sbmlstr)
flatcomb = antimony.getAntimonyString('combined')

r = te.loada(flatcomb)
r.exportToAntimony('combined.txt')
r.draw()