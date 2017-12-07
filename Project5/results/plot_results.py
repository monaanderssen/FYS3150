import numpy as np 
import matplotlib.pyplot as plt 
import sys

inputfile = sys.argv[1]
file = open(inputfile, "r")

numberOfIterations = []
time = []
T = []
E_k = []
E_p = []
E_tot = []
D = []
for line in file:
	line_ = line.split()
	numberOfIterations.append(float(line_[0]))
	time.append(float(line_[1]))
	T.append(float(line_[2]))
	E_k.append(float(line_[3]))
	E_p.append(float(line_[4]))
	E_tot.append(float(line_[5]))
	D.append(float(line_[6])) # we divide by time in the program and needs to include this factor. The other factor is for the length squared

plt.plot(time, T)
plt.show()

plt.plot(time, E_k, "b", time, E_p, "r", time, E_tot, "k")
plt.show()

plt.plot(time, D)
plt.title('Diffusion')
plt.xlabel('Time [seconds]')
plt.ylabel('Diffusion [m^2/s]')
plt.show()


