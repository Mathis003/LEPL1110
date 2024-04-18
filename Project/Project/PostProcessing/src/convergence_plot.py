# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import linregress

# # Charger les données
# data = np.genfromtxt('../../Processing/data/def_convergence.csv', delimiter=',', skip_header=1)

# # Extraire les données
# h = data[:, 0]
# nNodes = data[:, 1]
# error_X = data[:, 2]
# error_Y = data[:, 3]

# # Calculer l'erreur totale
# error_total = np.sqrt(error_X**2 + error_Y**2)

# # Définir l'erreur de référence
# ref_error_total = error_total[0]

# # Supprimer l'erreur de référence et la valeur correspondante de h
# error_total = error_total[1:]
# h = h[1:]

# # Calculer l'erreur absolue
# error_total = abs(ref_error_total - error_total)

# log_h = np.log(h)
# log_errors = np.log(error_total)
# slope, intercept, _, _, _ = linregress(log_h, log_errors)
# ordre_convergence = abs(slope)
# error_order = 1/ordre_convergence
# print(error_total)
# print(ordre_convergence)
# print(error_order)
# # Tracer le graphique
# plt.figure()

# # Plot with inverted x-axis
# plt.semilogx(h, error_total, marker='o', linestyle='-', color='b')
# plt.plot(h, np.exp(intercept) * h**slope, linestyle='--', color='r', label="Ordre de convergence: {:.2f}\nOrdre d'erreur: {:.2f}".format(ordre_convergence, error_order))
# plt.gca().invert_xaxis()

# # Labeling
# plt.xlabel('Taille du maillage (h)')
# plt.ylabel('Erreur de convergence')
# plt.title('Erreur de convergence en fonction du nombre de nœuds')
# plt.grid(True)
# plt.yscale("log")
# plt.savefig('../../Processing/data/convergence_plot.pdf', bbox_inches='tight')

# plt.legend()
# plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Charger les données
data = np.genfromtxt('../../Processing/data/def_convergence.csv', delimiter=',', skip_header=1)

# Extraire les données
h = data[:, 0]
nNodes = data[:, 1]
error_X = data[:, 2]
error_Y = data[:, 3]

# Calculer l'erreur totale
error_total = np.sqrt(error_X**2 + error_Y**2)

# Définir l'erreur de référence
ref_error_total = error_total[0]

# Supprimer l'erreur de référence et la valeur correspondante de h
error_total = error_total[1:]
h = h[1:]

# Calculer l'erreur absolue
error_total = abs(ref_error_total - error_total)

# Diviser les données en trois parties
h_part1, h_part2, h_part3 = h[:3], h[2:6], h[5:]
error_part1, error_part2, error_part3 = error_total[:3], error_total[2:6], error_total[5:]

# Calculer les régressions linéaires pour chaque partie
slope1, intercept1, _, _, _ = linregress(np.log(h_part1), np.log(error_part1))
slope2, intercept2, _, _, _ = linregress(np.log(h_part2), np.log(error_part2))
slope3, intercept3, _, _, _ = linregress(np.log(h_part3), np.log(error_part3))

# Calculer les ordres de convergence et d'erreur pour chaque partie
ordre_convergence1 = abs(slope1)
ordre_convergence2 = abs(slope2)
ordre_convergence3 = abs(slope3)
error_order1 = 1 / ordre_convergence1
error_order2 = 1 / ordre_convergence2
error_order3 = 1 / ordre_convergence3

plt.figure()

plt.semilogx(h, error_total, marker='o', linestyle='-', color='b')
plt.gca().invert_xaxis()

plt.xlabel('Mesh Size ($h$)')
plt.ylabel('Convergence Error')
plt.title('Convergence Error vs. Mesh Size')
plt.grid(True)
plt.yscale("log")
plt.savefig('../../Processing/data/convergence_plot.pdf', bbox_inches='tight')

# Use LaTeX formatting for legend
plt.legend(['Experimental Convergence Error'], loc='upper right', prop={'size': 10})
plt.show()