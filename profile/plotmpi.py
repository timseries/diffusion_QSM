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
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

def normvector(n): 
    nnew=[];
    chain=itertools.chain(*n)
    maxval=max(list(chain))
    for j in n:
        nnew=nnew+[[i/maxval for i in j]]
    return nnew

fnames=[]
for root, dirs, filenames in os.walk(indir):
    fnames.extend(filenames)
    break

mpifiles=[f for f in filenames if (f.startswith('mpi') and not(f.endswith('.viz')))]
ompfiles=[f for f in filenames if (f.startswith('pom') and not(f.endswith('.viz')))]
hpmfiles=[f for f in filenames if (f.startswith('hpm') and not(f.endswith('.viz')))]

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

iteration=0
rects=[]
width=.1
#normalize each of these datatypes
num_calls=normvector(num_calls);
avg_bytes=normvector(avg_bytes);
time=normvector(time);

ind = np.arange(len(routines[0]))
fig = plt.figure()
ax = fig.add_subplot(111)
for routine in routines: #get the normalized retangles for each group, iterates across files
    rects=rects+[ax.bar(ind,tuple(time[iteration]),width,color=cm.jet(1.*iteration/len(ind)))]
    ind=ind+width                                     
    iteration+=1
#build the data series
ax.set_title("Time spent in MPI Calls for experiment: " + indir.split('/')[-1])    
ax.set_ylabel('Normalized Time')
ax.set_xlabel('Process')
ax.set_xticklabels(tuple(processes))
ax.legend(tuple([i[0] for i in rects]),tuple(routines[0]))

for rect in rects:
    autolabel(rect)

plt.show()
print time
