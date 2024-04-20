import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('../../Rapport/data/csv/deformation.csv', delimiter=',', skip_header=1)

h = data[:, 0]
nNodes = data[:, 1]
error_X = data[:, 2]
error_Y = data[:, 3]

# Calculate the total error
error_total = np.sqrt(error_X * error_X + error_Y * error_Y)

# Take the reference error and remove it from the error_total
ref_error_total = error_total[0]
error_total = error_total[1:]
h = h[1:]

error_total = (ref_error_total - error_total) * (ref_error_total - error_total)

plt.figure()
plt.plot(h, error_total, marker='o', linestyle='-', color='b')
plt.gca().invert_xaxis()
plt.xlabel('Mesh Size ($h$)')
plt.ylabel('Convergence Error')
plt.title('Convergence Error vs. Mesh Size')
plt.grid(True)
plt.xscale("log")
plt.yscale("log")
# plt.savefig('../../Rapport/data/convergence_plot.pdf', bbox_inches='tight')
plt.legend(['Experimental Convergence Error'], loc='upper right', prop={'size': 10})
plt.show()