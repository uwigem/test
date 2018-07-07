# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 19:45:34 2018

@author: Yoshi
"""
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import timeit
import time
class HaltException(Exception): pass
np.set_printoptions(linewidth=160)

r = te.loada('''
             model testmodel()
             
             //// Reactions ////
             
             //X0 is "resource" [boundary species]
             //repressed by S2
             J0: $X0 -> S1 ; Vm1 * X0 / (km1 + X0 *(1 + S2/k10))
             
             J1: S1 -> S2 ; S1/(0.5 + S1);
             
             //a second input with "resource"
             J2: $X1 -> S2 ; .1*X1/(.5+X1);
             
             //uses cofactor A and returns B
             J3: S2 + A -> S4 + B ; Vm2 * A * S2 / (k3 + A + S2)
             
             //uses cofactor B and returns A
             J4: S1 + B -> S3 + A ; Vm3 * B * S1 / (k4 + B + S1)
             
             //repressed by S4
             J5: S2 -> S3 ; Vm4 * S2 / (km6 + S2 *(1 + S4/k11))
             
             J6: S4 -> S3 ; k7*S4
             
             //Y1 is "output" to environment [boundary species]
             J7: S3 -> $Y1 ; k8*S3   

             //second output (waste)
             J8: S4 -> $W ; k9*S4
             

             //// Values ////
            
            //species concentrations
             X0 = 1; X1 = 1; 
             S1 = 0; S2 = 0; S3 = 0; S4 = 0; 
             Y1 = 1; 
             W = 1
            
             //reaction rates
             km1 = 1; k2 = 1; k3 = 1; k4 = 1; 
             k5 = 1; km6 = 1; k7 = 1; k8 = 1;
             k9 = 1; k10 = 1; k11 = 1
             
             //Vm values             
             Vm1 = 2; Vm2 = .1; Vm3 = .1; Vm4 = .2;

             //rxn co-factors               
             A = 1; B = 1;
            
     end
''')

a = r.draw(layout='circo')

#%% ## Functions ##

def getIDs(model):
    ##collects all the data from the model, and then outputs a numpy vector. On first run, it will also output index vector.
    ids = np.array(r.model.getReactionIds()+r.model.getBoundarySpeciesIds()+r.model.getFloatingSpeciesIds())
    print ('ids are: ')
    count = 0
    for i in ids:
        print (str(count) +': ' + (i))
        count +=1
    print ('What is output? [Pick one using its index]')
    output = input()
    print ('What are your variables? (Put \',\' in between index) [Pick 3 to 6]')
    variables = input()
    origIds = ids
    print ('What will NOT be collected? (Put \',\' in between index)')
    nullOut = input()    
    selIds = tuple([ids[output]])
    for i in variables:
        selIds = selIds + tuple([ids[i]])
    indIds = list(variables)
    nullOutIds = list(nullOut)
    indIds.insert(0,output)
    return (selIds,indIds,origIds,nullOutIds)
    ##output: one tuple one list and a constant

def collectData(model,origIds,nullOutIds):
    ##collects all the data from the model, and then outputs a numpy vector. 
    combinedVals = np.concatenate((np.array(r.model.getReactionRates()),np.array(r.model.getBoundarySpeciesConcentrations()),np.array(r.model.getFloatingSpeciesConcentrations())))
    selVals = []
    dataIds = list()
    indDataIds = list()
    for i in range(len(origIds)):
        if i not in nullOutIds:
            selVals.append(combinedVals[i])
            dataIds.append(origIds[i])
            indDataIds.append(i)
    return selVals,dataIds,indDataIds
#   output: list

# checks for paraMatrix and calls createParaMat if not found
# will also check the working directory for paraMatrix.npy and load it if found
def checkParaMat(paraMatrix,simList,simCount,selIdsCount,selIds):
    if os.path.exists('./data/paraMatrix.npy') and os.path.exists('./data/selIdsCount.npy') and os.path.exists('./data/simCount.npy') and os.path.exists('./data/selIds.npy'):
        if np.isnan(paraMatrix).all():  #matrix is blank
            lastsimCount = np.load('./data/simCount.npy')
            lastselIdsCount = np.load('./data/selIdsCount.npy')
            lastIds = tuple(np.load('./data/selIds.npy'))
            if lastsimCount == simCount and lastselIdsCount == selIdsCount and lastIds == selIds:
                paraMatrix = np.load('./data/paraMatrix.npy')
                print "Parameters Matrix Loaded (this run is identical to previous run.)"
                return paraMatrix
        else:                            #matrix is not blank
            lastsimCount = np.load('./data/simCount.npy')
            lastselIdsCount = np.load('./data/selIdsCount.npy')
            lastIds = tuple(np.load('./data/selIds.npy'))
            if lastsimCount == simCount and lastselIdsCount == selIdsCount and lastIds == selIds:
                print "(Parameters matrix for this run is identical and already generated)"
                return paraMatrix
    print 'Parameters Matrix not saved'
    paraMatrix = createParaMat(paraMatrix,simList,selIdsCount)
    np.save('./data/paraMatrix.npy', paraMatrix)
    np.save('./data/simCount.npy',simCount)
    np.save('./data/selIdsCount.npy',selIdsCount)
    np.save('./data/selIds.npy',selIds)
    print "Parameters Matrix Created and saved"
    return paraMatrix

def createParaMat(paraMatrix,simList,selIdsCount):
    #create a NaN matrix to place parameters into.
    #columns represent parameters used for each simulation
    #rows represent boundary species
    count = -1
    paramRangeNum = 8  #turn into logspace
    print('Parameter range: ')
    count = 0
    for i in simList:
        print (str(paramRangeNum)+'**'+str(i))
        count +=1    
    #create a large matrix of all possible parameter combinations
    print('Creating ParaMat...')
    val0 = (paramRangeNum**ph for ph in simList) #creates generator going from 10**-simCount to 10**simCount
    #ph is a placeholder
    if selIdsCount>= 1:
        for x in simList:
            V0 = next(val0)    
            
            if selIdsCount>= 2:
                val1 = (paramRangeNum**ph for ph in simList)
                for y in simList:
                    V1 = next(val1)  
                    
                    if selIdsCount>= 3:
                        val2 = (paramRangeNum**ph for ph in simList)
                        for z in simList:
                            V2 = next(val2)
                            #print(x,y,z)  #debug
                            if selIdsCount>= 4:
                                val3 = (paramRangeNum**ph for ph in simList)
                                for a in simList:
                                    V3 = next(val3)
    
                                    if selIdsCount>= 5:
                                        val4 = (paramRangeNum**ph for ph in simList)
                                        for b in simList:
                                            V4 = next(val4)
                                            
                                            if selIdsCount== 6:
                                                val5 = (paramRangeNum**ph for ph in simList)
                                                for c in simList:
                                                    V5 = next(val5)     
                                                    count += 1
                                                    paraMatrix[:,count]=np.array([V0,V1,V2,V3,V4,V5])
                                            else:
                                                count += 1
                                                paraMatrix[:,count]=np.array([V0,V1,V2,V3,V4])
                                    else:
                                        count += 1
                                        paraMatrix[:,count]=np.array([V0,V1,V2,V3])
                            else:
                                count += 1
                                paraMatrix[:,count]=np.array([V0,V1,V2])
                    else:
                        count += 1
                        paraMatrix[:,count]=np.array([V0,V1])
            else:
                count += 1
                paraMatrix[:,count]=np.array([V0])
    
    #todo: probably want to use recursion and cartesian product in the future
    
    #cut out NaN values
    paraMatrix = paraMatrix[:,~np.all(np.isnan(paraMatrix), axis=0)]
    return paraMatrix

def dataSlice(data,criteria,range,arrange):
    data
    
    return slice
#%% ################################ Main method ###########################

selIds = ()

### load in model ###
r.steadyState() # run the simulation once with default values to generate a model object
selIds,indIds,origIds,nullOutIds = getIDs(r.model) # retrieve species IDs
selIdsCount = len(indIds)-1

### initializations ###
simCount = 10   # parameters will go from -simCount to simCount
stepSize = 1    # increment size of parameter values
simList = np.arange(-simCount,simCount+stepSize,stepSize).tolist()
paraMatrix = np.full((selIdsCount,5000000),np.nan,dtype=np.float64) 
    # 5 million simulations -> assume we won't need more than this.
    # todo: make this general

# generate paraMatrix if not done yet (should only run once)
# tries to regenerate one from saved file if existing.
# what if the user wants to generate a new one? -> This is governed by parameter changes
paraMatrix = checkParaMat(paraMatrix,simList,simCount,selIdsCount,selIds)

#%% Run simulations : generate totalData 

# set up a NaN matrix to collect simulation results
totalData = np.full((len(origIds)-len(nullOutIds),len(paraMatrix.T)+1),np.nan,dtype=np.float64)

print('\nRunning '+str(np.size(paraMatrix))+' simulations...')
if  np.size(paraMatrix)>1000000:
    print('You might want to step away for a coffee break.')
raw_input("Press Enter to continue ...")
print('Starting simulation generation...')
start = np.floor(timeit.default_timer())
for idx,sim in enumerate(paraMatrix.T):
    # todo: generalize the first element picked in paraMatrix with selIds
    r.reset()
    count = 0
    for i in selIds:
        if i !=selIds[0]:  #yes, we *want* to skip the output species (first species in selIds)
            r.setValue(i,np.float64(paraMatrix[count,idx]))
            count += 1
    try: # some set of variables may end up singular
        r.steadyState()
        if idx == 1:
            runData,dataIds,indDataIds = collectData(r.model,origIds,nullOutIds)
        else:
            runData,_,_ = collectData(r.model,origIds,nullOutIds)
        totalData[:,idx]=runData
    except RuntimeError as e:
        pass # if singular, just keep going, leave that column as a NaN
    #print(r.Y1) #debug
    if idx%500==0:
        print(str(idx) + ' out of ' + str(len(paraMatrix.T)))
    timeDuration = np.floor(timeit.default_timer()-start)
    if timeDuration-start==60:
        print('test')
        print('1 minute has passed...')
    if timeDuration-start==120:
        print('2 minutes has passed...')        
    if timeDuration-start==180:
        print('3 minutes has passed...')
    if timeDuration-start==240:
        print('4 minutes has passed...')                
    if timeDuration-start==300:
        print('5 minutes have passed...')
    if timeDuration-start==600:
        print('10 minutes have passed... Is your model ok?')
        raise SystemExit
    if timeDuration-start==900:
        print('15 minutes have passed... Hmmm....')
    if timeDuration-start==1200:
        print('20 minutes have passed... Get a better laptop lol')
#    print totalData[:,count]
print('simulation generation done.')
stop = timeit.default_timer()
print (str((stop-start)/60) + ' minutes runtime')

#%% arranging totalData to a dataframe
totalData = totalData[:,~np.all(np.isnan(totalData), axis=0)]

# Make a DataFrame with only dataIds species, transpose to make columns into species
data = pd.DataFrame(data=totalData.T,
                    columns=dataIds)

print('******* Creating a 3-D surface graph *********')

print ('Species able to be graphed are: ')
count = 0
for i in dataIds:
    print (str(count) +': ' + (i))
    count +=1
while True:
    try:
        print ('What species are you holding? (at least one)')
        out = input()
    except NameError:
        print('Only accepts indexes!')
        continue
    break
hold1 = dataIds[out]

print ('Do you want to hold something else? [1. Yes/2. No]')
out = input()
if out==1:
    print ('What species are you holding?')
    out = input()
    hold2 = dataIds[out]
    
print ('What species are x-axis and y-axis?')
out = input()
xaxis = dataIds[out[0]]
yaxis = dataIds[out[1]]

val = data.at[int(data.shape[0]/2),hold1]  #arbitrarily select the middle row to sample number of hold1
boole = data[hold1]==val
print('Setting '+hold1+' to '+str(val))    
if np.all(~boole)==True:    #All elements are False ; FIX
    print('WARNING: No points to plot')
slicedData = data[boole]

#slicedData = data #debug
if 'hold2' in locals():
    val = data.at[int(slicedData.shape[0]/2),hold2]     #arbitrarily select the middle row to sample number of hold2
    boole = slicedData[hold2]==val
    print('Setting '+hold2+' to '+str(val))    
    if np.all(~boole)==True:    #FIX  elements are False
        try:
            raise HaltException('WARNING: No points to plot')
        except HaltException as h:
            print(h)     
    slicedData = slicedData[boole]
#slicedData = slicedData[(slicedData[selIds[0]]>0.1)]

x = np.log10(slicedData[xaxis])
y = np.log10(slicedData[yaxis]) 
z = slicedData[selIds[0]]

fig = plt.figure()
ax = plt.gca(projection='3d')
#zi = griddata((X, Y), Z, (X[None,:], Y[:,None]), method='cubic')

#xi = np.linspace(-10, 10)
#yi = np.linspace(-10, 10)
#
#X, Y = np.meshgrid(xi, yi)
#Z = griddata(x, y, z, xi, yi)
#
#surf = ax.plot_surface(X,Y,Z,
#        linewidth=0,cmap='jet')

p = ax.scatter(x,y,z,cmap='gnuplot')
#for i in zip(x,y,z):
#    ax.plot(x(i,2), yval, zrow, color = 'black')
#plt.show()

## add colorbar
#fig.colorbar(p)

ax.set_xlabel(xaxis + ' (log)')
ax.set_ylabel(yaxis + ' (log)')
ax.set_zlabel(selIds[0])
plt.show()

#slicedData = dataSlice(data,criteria,paramRange,arrange)
