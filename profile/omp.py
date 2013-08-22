#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import profile_utils as pu
from matplotlib import cm
import pylab
import pdb
import re

def plotomp(indir,ompfiles):
    processes=[]
    time=[]
    allsecnames=[]
    allsecdata=[]
    allseclabels=[]
#read in all of the data and parse
    for fname in ompfiles:
        with open(indir+'/'+fname) as f:
            data = f.read()
            data = data.split('\n')
            processes = processes+[fname.split('_')[-1]]
            f.close()
        sinds=[]
        einds=[]
        seclabels=[]
        for i in range(0,len(data)):
            match=re.search('file:',data[i])
            if match:
                seclabels.append(data[i].strip())
                sinds.append(i+2)
            match=re.search('MeanThreadtime',data[i])
            if match: 
                einds.append(i+1)
                #pull out each of the sections
        eind=0
        secnames=[]
        secdata=[]
        for sind in sinds:
            secnames+=[[row.split(':')[0].strip() for row in data[sind:einds[eind]]]]
            secdata+=[map(float,[row.split(':')[1].split()[0] for row in data[sind:einds[eind]]])]
            eind+=1
            #for now this is just for the first hardware thread
#find the master metric_name  list used for matching all of the other parameters
            #combbine the sections for this file (process)
        allsecnames+=[secnames]
        allsecdata+=[secdata]
        allseclabels+=[seclabels]
#    pdb.set_trace()
    master_metric_names=secnames[0]
#    master_metric_names=sorted(list(set(list(itertools.chain(*metric_names)))),key=str.lower)
#    master_metric_names=sorted(list(set(list(itertools.chain(*metric_names)))),key=str.lower)
#    master_metric_names=master_metric_names[]
    iteration=0
    rects=[]
    metrics=[]
    width=.07

#do one plot for each instr_section
#    pdb.set_trace()
    for instr_section in range(0,len(allsecdata[0])):
        ind = np.arange(len(allsecdata))
        fig = plt.figure(instr_section)
        ax = fig.add_subplot(111)

        metrics=[pu.reorder_by_keys(master_metric_names,allsecnames[process][instr_section],allsecdata[process][instr_section]) for process in range(0, len(allsecdata))]
    #normalize each metric across all processes, but preserve percent overhead and load imbalance
        metric=[]
        rects=[]
        for metric_index in range(0,len(metrics[0])):
            if (master_metric_names[metric_index]=='PercentOverhead' \
            or master_metric_names[metric_index]=='LoadImbalance'):
                metric=[metrics[process][metric_index] for process in range(0,len(metrics))]
            else:
                metric=pu.normvector([[metrics[process][metric_index] for process in range(0,len(metrics))]])[0]
            rects=rects+[ax.bar(ind,tuple(metric),width,color=cm.jet(1.*metric_index/len(master_metric_names)))]
            ind=ind+width                                     
#store values back in big table for future use...
            for process in range(0,len(metrics)):
                metrics[process][metric_index]=metric[process]
        #plot...
        ax.set_title("Normalized OMP metrics for each process: " + indir.split('/')[-1])    
        ax.set_ylabel('Normalized Value or Percentage, where applicable')
        ax.set_xlabel('Process')
        ax.set_xticks(ind-5*width)
        ax.set_xticklabels(tuple(processes))
        ax.legend(tuple([i[0] for i in rects]),tuple(master_metric_names))
        for rect in rects:
            pu.autolabel(rect)      
#    pdb.set_trace()
    pylab.ion()    
    plt.show()
#print time
    