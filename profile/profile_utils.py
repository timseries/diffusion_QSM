import numpy as np
import itertools

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

def normvector2(n): 
    nnew=[];
    maxval=max(n)
    for j in n:
        nnew=nnew+[[i/maxval for i in j]]
    return nnew

def reorder_by_keys(masterkeys,inputkeys,inputvals):
    #inds=where each of the funtion calls appears in the master, sorted list
    inds=[masterkeys.index(f) for f in inputkeys]
    retval=np.ones(len(masterkeys))*-.01
    retval[inds]=inputvals
    #want to return values corresponding to these indices in the inputvals list
    return retval
