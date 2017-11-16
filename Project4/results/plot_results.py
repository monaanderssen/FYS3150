from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter



input = sys.argv[1:]

length = np.zeros(len(input))
i = 0
for input_ in input:
	file = open(input_, "r")
	for line in file: 
	    length[i] += 1
	print length[i]
	i += 1
print length

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

