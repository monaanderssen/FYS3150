import numpy as np 
import matplotlib.pyplot as plt 
import sys

inputfile = sys.argv[]
file = open(inputfile, "r")

t = []
T = []
E_k = []
E_p = []
E_tot = []
for line in file:
	line_ = line.split()
	t.append(float(line_[0]))
	T.append(float(line_[1]))
	E_k.append(float(line_[2]))
	E_p.append(float(line_[3]))
	E_tot.append(float(line_[4]))

plt.plot(t, T)
plt.show()


plt.plot(t, E_k, "b", t, E_p, "r", t, E_tot, "k")
plt.show()


