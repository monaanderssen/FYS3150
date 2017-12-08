import numpy as np 
import matplotlib.pyplot as plt 
import sys
import glob
import os

average_ratio = 0
counter = 0
legend_name_l = []
path = "/home/pederbh/UiO/FYS4150/FYS3150/Project5/results/ratio_seed"
for filename in sorted(glob.glob(os.path.join(path, '*.txt'))):
    with open(filename) as my_file:
        ratio = []
        time = []
        legend_name = filename.split('/')[-1]
        legend_name = legend_name.split('.')[0]
        legend_name =legend_name[7:]
        legend_name_l.append("$T_i$ = " + legend_name + "K")
        
        for line in my_file:
            line_ = line.split()
            time.append(float(line_[0]))
            ratio.append(float(line_[1]))

        if counter > 0:
            plt.plot(time, ratio)
            plt.hold('on')
            #plt.legend(legend_name)

        average_ratio += ratio[-1]
        counter += 1

average_ratio /= counter
print "Average ratio: {}".format(average_ratio)
print counter
plt.grid('on')
plt.title('Ratio $T/T_{initial}$')
plt.xlabel('Time [seconds]')
plt.ylabel('$T/T_i$')
#plt.legend(legend_name_l)
plt.show()


