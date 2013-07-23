#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import profile_utils as pu
from matplotlib import cm
import pdb;

def plothpm(indir,hpmfiles):
    metric_names=[]
    processes=[]
    metrics=[]
    time=[]
    for fname in hpmfiles:
        with open(indir+'/'+fname) as f:
            data = f.read()
            data = data.split('\n')
            data = [data[i] for i in range(26,32) + range(33,37)]
            metrics = metrics + [map(float,[row.split()[0] for row in data])]
            processes = processes+[fname.split('_')[-1]]
            #for now this is just for the first hardware thread
            metric_names = metric_names + [[" ".join([row.split()[i] for i in range(1,len(row.split()))]) for row in data]]
#find the master metric_name  list used for matching all of the other parameters
    master_metric_names=metric_names[0]
#    master_metric_names=sorted(list(set(list(itertools.chain(*metric_names)))),key=str.lower)
#    master_metric_names=sorted(list(set(list(itertools.chain(*metric_names)))),key=str.lower)
#    master_metric_names=master_metric_names[]
    iteration=0
    rects=[]
    width=.07
    metrics=[pu.reorder_by_keys(master_metric_names,metric_names[i],metrics[i]) for i in range(len(metrics))]
#normalize each of these datatypes and recombine them load misses and total loads
    for metric_index in range(0,len(metrics[0])):
        inds=range(0,len(metrics))
        metric_norm=pu.normvector([[metrics[i][metric_index] for i in inds]])
        for mset in range(0,len(metrics)):
            metrics[mset][metric_index]=metric_norm[0][mset]
            
#    pdb.set_trace()
    ind = np.arange(len(metric_names))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for metric_name in master_metric_names: #get the normalized retangles for each group, iterates across files
        rects=rects+[ax.bar(ind,tuple([x[iteration] for x in metrics]),width,color=cm.jet(1.*iteration/len(master_metric_names)))]
        ind=ind+width                                     
        iteration+=1
    #build the data series
    ax.set_title("Normalized HPM metrics for each process: " + indir.split('/')[-1])    
    #ax.set_ylabel('Normalized Time')
    ax.set_ylabel('Time (s)')
    ax.set_xlabel('Process')
    #ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(processes))
    ax.legend(tuple([i[0] for i in rects]),tuple(master_metric_names))
    for rect in rects:
        pu.autolabel(rect)
#    pylab.ion()    
    plt.show()
#print time









