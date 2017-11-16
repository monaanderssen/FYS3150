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
"""
#inputfiles = sys.argv[1:]
file = open(inputfiles[0], "r")
length = 0
for line in file: 
    length += 1
length *= len(inputfiles)"""

#M_C = np.linspace(0,length,length)

T = []
E_average  = []
E_variance = []
M_average = []
M_variance_T = []
M_abs = []

file_random = sys.argv[1]
file_ordered = sys.argv[2]
E_random_average = []
E_ordered_average = []
M_random_abs = []
M_ordered_abs = []

file_r = open(file_random, "r")

length = 0
for line in file_r: 
    length += 1

M_C = np.linspace(0,1000000,length)
file_r = open(file_random, "r")

for line in file_r: 
	line_ = line.split() 
	T.append(float(line_[0]))           
	E_random_average.append(float(line_[1]))
	M_random_abs.append(float(line_[5]))

file_o = open(file_ordered, "r")

for line in file_o: 
	line_ = line.split() 
	T.append(float(line_[0]))           
	E_ordered_average.append(float(line_[1]))
	M_ordered_abs.append(float(line_[5]))

#print len(E_random_average), len(M_C)

"""
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
"""
"""
plt.plot(M_C, E_random_average, 'b', M_C, E_ordered_average, 'g')
plt.title('Mean energy with T = 2.4', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('<E>', fontsize = 22)
plt.legend(['Random', 'Ordered'], fontsize=20)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.grid('on')
plt.axis([0,5000,-2,-0.5])
plt.show()
"""
"""
plt.plot(M_C, M_random_abs, 'b', M_C, M_ordered_abs, 'g')
plt.title('Magnetization with T = 2.4', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('<|M|>', fontsize = 22)
plt.legend(['Random', 'Ordered'], fontsize=20)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.grid('on')
plt.axis([0,5000,0,1])
plt.show()"""

#HISTOGRAM
plt.hist(E_random_average[1000:])
plt.title("Mean energy with T = 1.0", fontsize = 22)
plt.xlabel('Mean energy', fontsize = 22)
plt.ylabel('Probability', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.grid('on')
plt.show()
"""
fig, ax = plt.subplots()
plt.plot(T, E_variance)
plt.show()

plt.plot(M_C, E_average, '-', color='black')
plt.title('Mean energy', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('Mean energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.axis([0,10000,-3,-1])
plt.show()

plt.plot(M_C, M_abs)
plt.title('Mean absolutevalue magnetization', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
plt.ylabel('Mean energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.axis([0,10000,0,2])
plt.show()

"""
