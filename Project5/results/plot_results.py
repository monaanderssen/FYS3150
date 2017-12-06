import numpy as np 
import matplotlib.pyplot as plt 
import sys

inputfile = sys.argv[1]
file = open(inputfile, "r")

t = []
dt = []
T = []
E_k = []
E_p = []
E_tot = []
D = []
for line in file:
	line_ = line.split()
	t.append(float(line_[0]))
	dt.append(float(line_[1]))
	T.append(float(line_[2]))
	E_k.append(float(line_[3]))
	E_p.append(float(line_[4]))
	E_tot.append(float(line_[5]))
	D.append(float(line_[6]))

plt.plot(t, T)
plt.show()

plt.plot(t, E_k, "b", t, E_p, "r", t, E_tot, "k")
plt.show()

plt.plot(t, D)
plt.show()


