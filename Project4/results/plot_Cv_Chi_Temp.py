from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter
from numpy import argmax
from decimal import *


# L = 40, 60, 80, 100, T = [2.0, 2.6], dT = 0.01
inputfiles = sys.argv[1:]
print len(inputfiles)

number_of_spins_list = [20, 40, 60, 80, 100]

def write_from_file(file_, number_of_spins):
	T = []           
	E = []
	C_v = []
	#M_average = []
	Chi = []
	M_abs = []
	file = open(file_, "r")
	for line in file:
		line_ = line.split() 
		T.append(float(line_[0]))           
		E.append(float(line_[1])/float(number_of_spins**2))	# forgot to devide by nspins**2 in C++
		C_v.append(float(line_[2]))
		#M_average.append(float(line_[3]))
		Chi.append(float(line_[4]))
		M_abs.append(float(line_[5]))
	return T, E, C_v, Chi, M_abs



# Plotting C_v as function of T
for i in range(len(inputfiles)):
	T, E, C_v, Chi, M_abs = write_from_file(inputfiles[i], number_of_spins_list[i])

	plt.plot(T, Chi)
	plt.hold('on')
plt.xlabel("Temperature", fontsize=22); plt.ylabel("Susceptibility", fontsize=22);
plt.title("Susceptibility with lattice size LxL", fontsize=22)
plt.legend(["20x20", "40x40", "60x60", "80x80", "100x100"], fontsize=18)
pylab.xticks(fontsize=16); pylab.yticks(fontsize=16)
plt.xlim(2, 2.6)
plt.grid('on')
plt.show()

# Plotting Chi as function of T
for i in range(len(inputfiles)):
	T, E, C_v, Chi, M_abs = write_from_file(inputfiles[i], number_of_spins_list[i])
	
	#Finding the temperature to corresponding max value of C_V
	max_C_v = max(C_v)  # Find the maximum y value
	max_T = T[np.array(C_v).argmax()]  # Find the x value corresponding to the maximum y value
	print Decimal(max_T), max_C_v
	
	plt.plot(T, C_v)
	plt.hold('on')
plt.xlabel("Temperature", fontsize=22); plt.ylabel("Specific heat capacity", fontsize=22);
plt.title("Phase transitions with lattice size LxL", fontsize=22)
plt.legend(["20x20", "40x40", "60x60", "80x80", "100x100"], fontsize=18)
pylab.xticks(fontsize=16); pylab.yticks(fontsize=16)
plt.xlim(2, 2.6)
plt.grid('on')
plt.show()


# Plotting E as function of T
for i in range(len(inputfiles)):
	T, E, C_v, Chi, M_abs = write_from_file(inputfiles[i], number_of_spins_list[i])
	
	plt.plot(T, E)
	plt.hold('on')
plt.xlabel("Temperature", fontsize=22); plt.ylabel("< E >", fontsize=22);
plt.title("Mean energy with lattice size LxL", fontsize=22)
plt.legend(["20x20", "40x40", "60x60", "80x80", "100x100"], fontsize=18, loc=2)
pylab.xticks(fontsize=16); pylab.yticks(fontsize=16)
plt.xlim(2, 2.6)
plt.show()


# Plotting M_abs as function of T
for i in range(len(inputfiles)):
	T, E, C_v, Chi, M_abs = write_from_file(inputfiles[i], number_of_spins_list[i])
	
	plt.plot(T, M_abs)
	plt.hold('on')
plt.xlabel("Temperature", fontsize=22); plt.ylabel("<|M|>", fontsize=22);
plt.title("Mean absolute magnetization with lattice size LxL", fontsize=22)
plt.legend(["20x20", "40x40", "60x60", "80x80", "100x100"], fontsize=18)
pylab.xticks(fontsize=16); pylab.yticks(fontsize=16)
plt.xlim(2, 2.6)
plt.show()














