#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

with open("mpi_profile_0_0") as f:
    data = f.read()

data = data.split('\n')

# pull out the summary statistics table
functions = [row.split(' ')[0] for row in data[4:15]]
num_calls = [row.split(' ')[1] for row in data[4:15]]
avg_bytes = [row.split(' ')[2] for row in data[4:15]]
time = [row.split(' ')[2] for row in data[4:15]]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("HPM")    
ax1.set_xlabel('your x label..')
ax1.set_ylabel('your y label...')

ax1.plot(x,y, c='r', label='the data')

leg = ax1.legend()

plt.show()

