# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018

@author: Yoshi
"""

import os
import glob
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np
import antimony
import string
import itertools
import re

np.set_printoptions(linewidth=160)

def product_gen(n):
    for r in itertools.count(1):
        for i in itertools.product(n, repeat=r):
            yield "".join(i)

#%% Test model: putting two models together
numReads = 2
reads = {}
#for i in np.arange(0,numReads):
reads[0] = te.readFromFile('FFL1_ant.txt')
reads[1] = te.readFromFile('FFL2_ant.txt')
    #generalize later
    #combine this with for loop below this later
# Note: The string object name has to match the model name

readsalpha = zip(string.ascii_uppercase[:numReads],reads.values())

models = reads[0]
species = ""
currStr = readsalpha[0][1]
splitStr = re.split('(model)( [\*]?)(.*?)\\n\\n',currStr)
modules = readsalpha[0][0] + " : " + splitStr[3] + "; \n"
inputs = ""
outputs = ""
for i in np.arange(1,numReads):
    models = models + reads[i]
    currStr = readsalpha[i][1]
    splitStr = re.split('(model)( [\*]?)(.*?)\\n\\n',currStr)
    modules = modules + readsalpha[i][0] + " : " + splitStr[3] + "; \n" #i
    
    species = species + "var species p_c"+str(i)+"; \n" #i-1
    
    inputs = inputs + readsalpha[i-1][0] + ".p_input" + " is p_c" + str(i) + "; \n" #i-1
    outputs = outputs + readsalpha[i][0] + ".p_output" + " is p_c" + str(i) + "; \n" #i-1

species.replace('\\n', '\n')
modules.replace('\\n', '\n')
inputs.replace('\\n', '\n')
outputs.replace('\\n', '\n')

combined = models + 'model combined '+ modules + species + inputs + outputs + 'end'
    
##What I want
#combined = add(all models I have) + '''
#model combined
#    #repeat for n-1 modules
#    #!!module names have to be identical to model names from imported models!!
#    A : model1();
#    B : model2();
#    #repeat for all models I have
#    #specify a global input node
#    var species p_c;
#    A.p_input is p_c;
#    B.p_output is p_c;
#    #repeat for n-1 modules
# end 
# '''

#properly flatten the combined model
antimony.loadAntimonyString(combined)
sbmlstr = antimony.getSBMLString('combined')
antimony.loadSBMLString(sbmlstr)
flatcomb = antimony.getAntimonyString('combined')

r = te.loada(flatcomb)
r.exportToAntimony('combined.txt')
r.draw(layout='dot')


#combined =  reads +  '''
#model combined
#    var species p_c;
#    #!!module names have to be identical to model names from imported models!!
#    A : ffl1();
#    B : ffl1();
#    A.p_input is p_c;
#    B.p_output is p_c;
#end
#'''