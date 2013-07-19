#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib import cm

indir='../profile_results/07182013/n1p8t8openmp'

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        #ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
         #       ha='center', va='bottom')

def normvector(n): 
    nnew=[];
    chain=itertools.chain(*n)
    maxval=max(list(chain))
    for j in n:
        nnew=nnew+[[i/maxval for i in j]]
    return nnew

def reorder(masterkeys,inputkeys,inputvals):
    #inds=where each of the funtion calls appears in the master, sorted list
    inds=[masterkeys.index(f) for f in inputkeys]
    retval=np.ones(len(masterkeys))*-.01
    retval[inds]=inputvals
    #want to return values corresponding to these indices in the inputvals list
    return retval

fnames=[]
for root, dirs, filenames in os.walk(indir):
    fnames.extend(filenames)
    break

mpifiles=sorted([f for f in filenames if (f.startswith('mpi') and not(f.endswith('.viz')))],
                key=str.lower)
ompfiles=sorted([f for f in filenames if (f.startswith('pom') and not(f.endswith('.viz')))],
                key=str.lower)
hpmfiles=sorted([f for f in filenames if (f.startswith('hpm') and not(f.endswith('.viz')))],
                key=str.lower)

iteration=1    
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
num_calls=normvector(num_calls)
avg_bytes=normvector(avg_bytes)
time=normvector(time)

#reorder these to match the master_routines
num_calls=[reorder(master_routines,routines[i],num_calls[i]) for i in range(len(num_calls))]
avg_bytes=[reorder(master_routines,routines[i],avg_bytes[i]) for i in range(len(avg_bytes))]
time=[reorder(master_routines,routines[i],time[i]) for i in range(len(time))]

ind = np.arange(len(routines))
fig = plt.figure()
ax = fig.add_subplot(111)
for routine in master_routines: #get the normalized retangles for each group, iterates across files
    rects=rects+[ax.bar(ind,tuple([x[iteration] for x in time]),width,color=cm.jet(1.*iteration/len(master_routines)))]
    ind=ind+width                                     
    iteration+=1
#build the data series
ax.set_title("Time spent in MPI Calls for experiment: " + indir.split('/')[-1])    
ax.set_ylabel('Normalized Time')
ax.set_xlabel('Process')
#ax.set_xticks(ind+width)
ax.set_xticklabels(tuple(processes))
ax.legend(tuple([i[0] for i in rects]),tuple(master_routines))

for rect in rects:
    autolabel(rect)

plt.show()
#print time
