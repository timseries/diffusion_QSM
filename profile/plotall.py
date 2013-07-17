#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt

indir='../profile_results/07182013/n1p8t8openmp'


for root, dirs, filenames in os.walk(indir):
    mpifiles=[f for f in filenames  if (f.startswith('mpi') and not(f.endswith('.viz')))]
    ompfiles=[f for f in filenames  if (f.startswith('pom') and not(f.endswith('.viz')))]
    hpmfiles=[f for f in filenames  if (f.startswith('hpm') and not(f.endswith('.viz')))]
    
for f in mpifiles:


# with open("mpi_profile_0_0") as f:
#     data = f.read()

# data = data.split('\n')

# # pull out the summary statistics table
# labels = [row.split()[0] for row in data[4:15]]
# num_calls = map(float,[row.split()[1] for row in data[4:15]])
# avg_bytes = map(float,[row.split()[2] for row in data[4:15]])
# time = map(float,[row.split()[3] for row in data[4:15]])

# fig = plt.figure()

# ax1 = fig.add_subplot(111)
# ind = np.arrange(len(labels))
# rects1 = ax.bar(ind,
# ax1.set_title("Plot title...")    
# ax1.set_xlabel('your x label..')
# ax1.set_ylabel('your y label...')

# ax1.plot(x,y, c='r', label='the data')

# leg = ax1.legend()

# plt.show()

