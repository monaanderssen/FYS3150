from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter


# L = 40, 60, 80, 100, T = [2.0, 2.6], dT = 0.01
inputfiles = sys.argv[1:]
print len(inputfiles)


def write_from_file(file_):
	T = []           
	#E_average = []
	C_v = []
	#M_average = []
	Chi = []
	#M_abs = []

	file = open(file_, "r")
	for line in file:
		line_ = line.split() 
		T.append(float(line_[0]))           
		#E_average.append(float(line_[1]))
		C_v.append(float(line_[2]))
		#M_average.append(float(line_[3]))
		Chi.append(float(line_[4]))
		#M_abs.append(float(line_[5]))
	return T, C_v, Chi

# Plotting C_v as function of T
for i in range(len(inputfiles)):
	T, C_v, Chi = write_from_file(inputfiles[i])

	plt.plot(T, Chi)
	plt.hold('on')
plt.show()

# Plotting Chi as function of T
for i in range(len(inputfiles)):
	T, C_v, Chi = write_from_file(inputfiles[i])

	plt.plot(T, C_v)
	plt.hold('on')
plt.show()
