#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import profile_utils as pu
from matplotlib import cm

def plotmpi(indir,mpifiles):
    routines=[]
    processes=[]
    num_calls=[]
    avg_bytes=[]
    time=[]
    for fname in mpifiles:
        with open(indir+'/'+fname) as f:
            data = f.read()
            data = data.split('\n')
            routines = routines+[[row.split()[0] for row in data[4:15]]]
            processes = processes+[fname.split('_')[-1]]
            num_calls = num_calls+[map(float,[row.split()[1] for row in data[4:15]])]
            avg_bytes = avg_bytes+[map(float,[row.split()[2] for row in data[4:15]])]
            time = time+[map(float,[row.split()[3] for row in data[4:15]])]
#find the master routine list used for matching all of the other parameters
        
    master_routines=sorted(list(set(list(itertools.chain(*routines)))),key=str.lower)
    iteration=0
    rects=[]
    width=.07
#normalize each of these datatypes
#num_calls=normvector(num_calls)
#avg_bytes=normvector(avg_bytes)
#time=normvector(time)
#reorder these to match the master_routines
    num_calls=[pu.reorder_by_keys(master_routines,routines[i],num_calls[i]) for i in range(len(num_calls))]
    avg_bytes=[pu.reorder_by_keys(master_routines,routines[i],avg_bytes[i]) for i in range(len(avg_bytes))]
    time=[pu.reorder_by_keys(master_routines,routines[i],time[i]) for i in range(len(time))]
    ind = np.arange(len(routines))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for routine in master_routines: #get the normalized retangles for each group, iterates across files
        rects=rects+[ax.bar(ind,tuple([x[iteration] for x in time]),width,color=cm.jet(1.*iteration/len(master_routines)))]
        ind=ind+width                                     
        iteration+=1
    #build the data series
    ax.set_title("Time spent in MPI Calls for experiment: " + indir.split('/')[-1])    
    #ax.set_ylabel('Normalized Time')
    ax.set_ylabel('Time (s)')
    ax.set_xlabel('Process')
    #ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(processes))
    ax.legend(tuple([i[0] for i in rects]),tuple(master_routines))
    for rect in rects:
        pu.autolabel(rect)
#    pylab.ion()    
    plt.show(block=True)
#print time
