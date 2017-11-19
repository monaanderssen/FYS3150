from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter



inputfiles = sys.argv[1:]

def ConfigurationsAndCycles(inputfiles):
	length = np.zeros(len(inputfiles))
	i = 0
	for input_ in inputfiles:
		file = open(input_, "r")
		for line in file: 
		    length[i] += 1
		print length[i]
		i += 1
	
	def configurations_(length, filename):
		cycles = np.zeros(length)
		configurations = np.zeros(length)
		
		file = open(filename, "r")
		counter = 0
		for line in file: 
		    line_ = line.split()
		    cycles[counter] = float(line_[0])                
		    configurations[counter] = float(line_[1])
	
		    counter += 1
		return cycles, configurations
	
	for i in range(len(input)):
		cycles, configurations = configurations_(length[i], input[i])
	
		plt.plot(cycles, configurations, '-')
		plt.hold('on')
	
	plt.legend(['T = 1.0','T = 1.1','T = 1.2'])
	plt.title('Accepted configurations', fontsize = 22)
	plt.xlabel('Monte Carlo simulations', fontsize = 22)
	#ax.ticklabel_format(useOffset=False)
	plt.ylabel('Configurations', fontsize = 22)# change to v_tilde(x)
	#pylab.xticks(fontsize=16)
	#pylab.yticks(fontsize=16)
	#plt.legend(['Sun'])
	plt.show()

#ConfigurationAndCycles(input)

def Chi_temp(inputfiles):
	T = []           
	E_average = []
	E_variance = []
	M_average = []
	M_variance_T = []
	M_abs = []
	length = 0
	for input_ in inputfiles:
		file = open(input_, "r")
		for line in file: 
		    length += 1
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
	plt.plot(T, M_variance_T)
	plt.show()



Chi_temp(inputfiles)





