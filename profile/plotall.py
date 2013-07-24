#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import mpi
import hpm
import omp
import profile_utils as pu
import pylab
from matplotlib import cm
import getopt
import sys
def main():
    opts, extraparams = getopt.getopt(sys.argv[1:],'d:',["directory"])
    print opts
    for o,p in opts:
        if o in ('-d','--directory'):
            indir=p
    #indir='../profile_results/07182013/n1p8t8openmp'
    fnames=[]
    print opts
    print indir
    for root, dirs, filenames in os.walk(indir):
        fnames.extend(filenames)
        break
    mpifiles=sorted([f for f in filenames if (f.startswith('mpi') and not(f.endswith('.viz')))],
                key=str.lower)
    ompfiles=sorted([f for f in filenames if (f.startswith('pom') and not(f.endswith('.viz')))],
                key=str.lower)
    hpmfiles=sorted([f for f in filenames if (f.startswith('hpm') and not(f.endswith('.viz')))],
                key=str.lower)
    pylab.ion()    
#    omp.plotomp(indir,ompfiles)
 #   hpm.plothpm(indir,hpmfiles)
    mpi.plotmpi(indir,mpifiles)
    return 1
if __name__ == "__main__":
    main()
