import numpy as np
import matplotlib.pyplot as plt

precision, nodes, full, no_order, x_order, y_order, rmck_order = np.loadtxt('../../Processing/data/solver_type.csv', delimiter=',', unpack=True, skiprows=1)

plt.figure("Performance des différents solveurs")
plt.plot(nodes, full, '--o', label='Full')
plt.plot(nodes, no_order, '--o', label='No order')
plt.plot(nodes, x_order, '--o', label='X order')
plt.plot(nodes, y_order, '--o', label='Y order')
plt.plot(nodes, rmck_order, '--o', label='RMCK order')
plt.xlabel('Nombre de noeuds')
plt.ylabel('Temps (s)')
plt.legend()
plt.yscale('log')
plt.show()