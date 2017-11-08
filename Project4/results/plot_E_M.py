from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter

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
    M_abs[counter] = float(line_[3])

    counter += 1

fig, ax = plt.subplots()

plt.plot(M_C, E_average, '-', color='black')

plt.title('Mean energy', fontsize = 22)
plt.xlabel('Monte Carlo simulations', fontsize = 22)
#ax.ticklabel_format(useOffset=False)
plt.ylabel('Mean energy', fontsize = 22)
#pylab.xticks(fontsize=16)
#pylab.yticks(fontsize=16)
#plt.legend(['Sun'])
plt.show()


