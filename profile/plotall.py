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

def main():
    indir='../profile_results/07182013/n1p8t8openmp'
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
    pylab.ion()    
    omp.plotomp(indir,ompfiles)
    hpm.plothpm(indir,hpmfiles)
    mpi.plotmpi(indir,mpifiles)
    return 1
if __name__ == "__main__":
    main()
