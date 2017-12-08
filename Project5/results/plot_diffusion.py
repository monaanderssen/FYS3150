import numpy as np 
import matplotlib.pyplot as plt 
import sys

inputfile = sys.argv[1]
file = open(inputfile, "r")

D = []
T = []

for line in file:
    line_ = line.split()
    T.append(float(line_[0]))
    D.append(float(line_[1]))

plt.plot(T, D, 'o')
plt.grid('on')
plt.title('Diffusion')
plt.xlabel('Equilibrium temperature [K]')
plt.ylabel('Diffusion [m^2/s]')
plt.show()


