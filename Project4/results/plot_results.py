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

print length
T = np.zeros(length)
E_average = np.zeros(length)
C_V = np.zeros(length)
M = np.zeros(length)
chi = np.zeros(length)
M_abs = np.zeros(length)

file = open(input, "r")
counter = 0
for line in file: 
    line_ = line.split()
    T[counter] = float(line_[0])                
    E_average[counter] = float(line_[1])
    C_V[counter] = float(line_[2])
    M[counter] = float(line_[3])
    chi[counter] = float(line_[4])
    M_abs[counter] = float(line_[5])

    counter += 1

print C_V

fig, ax = plt.subplots()

plt.plot(T[:], C_V[:], '-', color='black')

plt.title('Heat capacity', fontsize = 22)
plt.xlabel('T', fontsize = 22)
#ax.ticklabel_format(useOffset=False)
plt.ylabel('C_V', fontsize = 22)# change to v_tilde(x)
#pylab.xticks(fontsize=16)
#pylab.yticks(fontsize=16)
#plt.legend(['Sun'])
plt.show()


