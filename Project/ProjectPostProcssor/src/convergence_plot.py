import numpy as np
import matplotlib.pyplot as plt

# !! LOST THE DATAS !!
# Forgot to push the csv file to the repo, and then "git stash", rest in peace "deformation_beam.csv" :(
# This data is false because we did some error in the calculation of the error in the code
# Check the real graph in the rapport (with the real data)

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

# error_total = (ref_error_total - error_total) * (ref_error_total - error_total)
error_total = (ref_error_total - error_total) * (ref_error_total - error_total)

plt.figure()
plt.plot(h, error_total, marker='o', linestyle='-', color='b', label='Experimental Convergence Error')
plt.plot(h, (h**2 / 50000)**2, linestyle='--', color='r', label='Convergence Order: 2')
plt.plot(h, (h**0.75 / 80000)**2, linestyle='--', color='g', label='Convergence Order: 0.75')
plt.xlabel('Mesh Size ($h$)')
plt.ylabel('Convergence Error')
plt.title('Convergence Error vs. Mesh Size')
plt.grid(True)
plt.xscale("log")
plt.yscale("log")
plt.legend(loc='upper left', prop={'size': 10})
# plt.savefig('../../Rapport/data/convergence_plot.pdf', bbox_inches='tight')
plt.show()