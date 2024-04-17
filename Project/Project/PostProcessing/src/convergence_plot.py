import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

data = np.genfromtxt('../../Processing/data/def_convergence.csv', delimiter=',', skip_header=1)

h = data[:, 0]
nNodes = data[:, 1]
error_X = data[:, 2]
error_Y = data[:, 3]

error_total = np.sqrt(error_X**2 + error_Y**2)

ref_error_total = error_total[0]
error_total = error_total[1:]
h = h[1:]
# error_total = (ref_error_total - error_total) * (ref_error_total - error_total)

error_total = abs(ref_error_total - error_total)

# 3.84, 256, -0.4519315326577248437267542158224387094378471374511718750000000000, -0.7117170313810334825177505990723147988319396972656250000000000000

slope, intercept, r_value, p_value, std_err = linregress(np.log(h), np.log(error_total))
order_error = slope
order_convergence = -1 / slope

print("Ordre d'erreur:", order_error)
print("Ordre de convergence:", order_convergence)

plt.figure()

plt.plot(h, np.exp(intercept) * h**slope, linestyle='--', color='r', label=f'Régression linéaire (pente={slope:.2f})')
plt.plot(h, error_total, marker='o', linestyle='-', color='b')
plt.xlabel('Taille du maillage (h)')
plt.ylabel('Erreur de convergence')
plt.title('Erreur de convergence en fonction du nombre de nœuds')
plt.grid(True)
plt.yscale('log')
plt.xscale('log')
plt.savefig('../../Processing/data/convergence_plot.pdf')
plt.legend()
plt.show()
