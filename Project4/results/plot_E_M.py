from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter
"""
input = sys.argv[1]

file = open(input, "r")

length = 0
for line in file: 
    length += 1

M_C = np.linspace(10,10000, length)
E_average = np.zeros(length)
M_abs = np.zeros(length)

file = open(input, "r")
counter = 0
for line in file: 
    line_ = line.split()            
    E_average[counter] = float(line_[1])
    M_abs[counter] = float(line_[5])

    counter += 1
"""
inputfiles = sys.argv[1:]
file = open(inputfiles[0], "r")
length = 0
for line in file: 
    length += 1
length *= len(inputfiles)

M_C = np.linspace(10,10000, length)
T = []
E_average  = []
E_variance = []
M_average = []
M_variance_T = []
M_abs = []

for files in inputfiles:
	file = open(files, "r")
	for line in file:
		line_ = line.split() 
		T.append(float(line_[0]))           
		E_average.append(float(line_[1]))
		E_variance.append(float(line_[2]))
		M_average.append(float(line_[3]))
		M_variance_T.append(float(line_[4]))
		M_abs.append(float(line_[5]))


fig, ax = plt.subplots()

plt.plot(M_C, E_average, '-', color='black')
plt.title('Mean energy', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('Mean energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.show()

plt.plot(M_C, M_abs)
plt.title('Mean absolutevalue magnetization', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('Mean energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.show()


