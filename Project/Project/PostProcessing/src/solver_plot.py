# import numpy as np
# import matplotlib.pyplot as plt

# mesh_size, nodes, full_solver, band_solver_no_order, band_solver_x_order, band_solver_y_order, band_solver_rcmk_order, no_order_ratio, x_order_ratio, y_order_ratio, rcmk_order_ratio = np.loadtxt('../../Processing/data/solver_type.csv', delimiter=',', unpack=True, skiprows=1)

# size = 2 * nodes
# band_no_order   = size / no_order_ratio
# band_x_order    = size / x_order_ratio
# band_y_order    = size / y_order_ratio
# band_rcmk_order = size / rcmk_order_ratio

# average_band_ratios = [np.mean(no_order_ratio), np.mean(x_order_ratio), np.mean(y_order_ratio), np.mean(rcmk_order_ratio)]
# reorder_labels      = ['No reordering', 'X reordering', 'Y reordering', 'RCMK reordering']

# complexity_full_solver   = np.power(nodes, 3) / 100000000
# complexity_band_solver   = np.power(nodes, 1.65) / 5000000
# complexity_no_reordering = size / 1.3
# complexity_x_reordering  = nodes ** 0.5 * 13

# plt.figure("Performance of Different Solvers with Different Reordering")
# plt.title(r"Performance of Different Solvers with Different Reordering")
# plt.plot(nodes, full_solver, '--o', label='Full Solver')
# plt.plot(nodes, band_solver_no_order, '--o', label='Band Solver without reordering')
# plt.plot(nodes, band_solver_x_order, '--o', label='Band Solver with X reordering')
# plt.plot(nodes, band_solver_y_order, '--o', label='Band Solver with Y reordering')
# plt.plot(nodes, band_solver_rcmk_order, '--o', label='Band Solver with RCMK reordering')
# plt.plot(nodes, complexity_full_solver, '-', label=r'Complexity of Full Solver: $O(n^3)$')
# plt.plot(nodes, complexity_band_solver, '-', label=r'Complexity of Band Solver: $O(n^{1.65})$')
# plt.xlabel(r'Number of Nodes [$n$]')
# plt.ylabel(r'Time [s]')
# plt.yscale('log')
# plt.xscale('log')
# plt.legend(loc='upper left', fontsize='small')


# # Adding average band ratios to the plot as a box in the bottom right corner
# # bbox_props = dict(boxstyle="round", facecolor="white", alpha=0.5)
# # plt.text(0.95, 0.05, '\n'.join([f'{label}: {ratio:.2f}' for label, ratio in zip(reorder_labels, average_band_ratios)]),
# #          transform=plt.gca().transAxes, verticalalignment='bottom', horizontalalignment='right', bbox=bbox_props)
# # plt.text(0.95, 0.24, 'Average Band Ratios :', transform=plt.gca().transAxes, verticalalignment='bottom', horizontalalignment='right')

# plt.savefig('../../Processing/data/solver_plot.svg')
# plt.savefig('../../Processing/data/solver_plot.pdf')

# plt.show()


# plt.figure("Band Size with Different Reordering")
# plt.title(r"Band Size with Different Reordering")
# plt.plot(nodes, band_no_order, '--o', label='Band Size without reordering')
# plt.plot(nodes, band_x_order, '--o', label='Band Size with X reordering')
# plt.plot(nodes, band_y_order, '--o', label='Band Size with Y reordering')
# plt.plot(nodes, band_rcmk_order, '--o', label='Band Size with RCMK reordering')
# plt.plot(nodes, complexity_no_reordering, '-', label=r'Complexity of Band Size without reordering: $O(n)$')
# plt.plot(nodes, complexity_x_reordering, '-', label=r'Complexity of Band Size with reordering: $O(\sqrt{n})$')
# plt.xlabel(r'Number of Nodes [$n$]')
# plt.ylabel(r'Band size [$\#$]')
# plt.yscale('log')
# plt.xscale('log')
# plt.legend(loc='upper left', fontsize='small')

# plt.savefig('../../Processing/data/band_size_plot.svg')
# plt.savefig('../../Processing/data/band_size_plot.pdf')

# plt.show()


import numpy as np
import matplotlib.pyplot as plt

# Chargement des données
mesh_size, nodes, full_solver, band_solver_no_order, band_solver_x_order, band_solver_y_order, band_solver_rcmk_order, no_order_ratio, x_order_ratio, y_order_ratio, rcmk_order_ratio = np.loadtxt('../../Processing/data/solver_type.csv', delimiter=',', unpack=True, skiprows=1)

# Calcul des valeurs nécessaires
size = 2 * nodes
band_no_order   = size / no_order_ratio
band_x_order    = size / x_order_ratio
band_y_order    = size / y_order_ratio
band_rcmk_order = size / rcmk_order_ratio

average_band_ratios = [np.mean(no_order_ratio), np.mean(x_order_ratio), np.mean(y_order_ratio), np.mean(rcmk_order_ratio)]
reorder_labels      = ['No reordering', 'X reordering', 'Y reordering', 'RCMK reordering']

complexity_full_solver   = np.power(nodes, 3) / 100000000
complexity_band_solver   = np.power(nodes, 1.65) / 5000000
complexity_no_reordering = size / 1.3
complexity_x_reordering  = nodes ** 0.5 * 13

# Création des sous-graphiques
fig, axes = plt.subplots(2, 1, figsize=(8, 10))

# Premier sous-graphique : Performance des solveurs avec différentes renumérotations
axes[0].set_title("Performance of Different Solvers with Different Reordering")
axes[0].plot(nodes, full_solver, '--o', label='Full Solver')
axes[0].plot(nodes, band_solver_no_order, '--o', label='Band Solver without reordering')
axes[0].plot(nodes, band_solver_x_order, '--o', label='Band Solver with X reordering')
axes[0].plot(nodes, band_solver_y_order, '--o', label='Band Solver with Y reordering')
axes[0].plot(nodes, band_solver_rcmk_order, '--o', label='Band Solver with RCMK reordering')
axes[0].plot(nodes, complexity_full_solver, '-', label=r'Complexity of Full Solver: $O(n^3)$')
axes[0].plot(nodes, complexity_band_solver, '-', label=r'Complexity of Band Solver: $O(n^{1.65})$')
axes[0].set_xlabel(r'Number of Nodes [$n$]')
axes[0].set_ylabel(r'Time [s]')
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].legend(loc='upper left', fontsize='small')

# Deuxième sous-graphique : Taille de la bande avec différentes renumérotations
axes[1].set_title("Band Size with Different Reordering")
axes[1].plot(nodes, band_no_order, '--o', label='Band Size without reordering')
axes[1].plot(nodes, band_x_order, '--o', label='Band Size with X reordering')
axes[1].plot(nodes, band_y_order, '--o', label='Band Size with Y reordering')
axes[1].plot(nodes, band_rcmk_order, '--o', label='Band Size with RCMK reordering')
axes[1].plot(nodes, complexity_no_reordering, '-', label=r'Complexity of Band Size without reordering: $O(n)$')
axes[1].plot(nodes, complexity_x_reordering, '-', label=r'Complexity of Band Size with reordering: $O(\sqrt{n})$')
axes[1].set_xlabel(r'Number of Nodes [$n$]')
axes[1].set_ylabel(r'Band size [$\#$]')
axes[1].set_yscale('log')
axes[1].set_xscale('log')
axes[1].legend(loc='upper left', fontsize='small')

plt.tight_layout()  # Pour éviter que les labels se chevauchent
plt.savefig('../../Processing/data/solver_plot.pdf')
plt.show()
