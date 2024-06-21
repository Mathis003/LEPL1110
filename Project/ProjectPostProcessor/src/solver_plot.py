import numpy as np
import matplotlib.pyplot as plt

mesh_size, nodes, full_solver, band_solver_no_order, band_solver_x_order, band_solver_y_order, band_solver_rcmk_order, no_order_ratio, x_order_ratio, y_order_ratio, rcmk_order_ratio = np.loadtxt('../../Rapport/data/csv/solver_type.csv', delimiter=',', unpack=True, skiprows=1)

size = 2 * nodes
band_no_order   = size / no_order_ratio
band_x_order    = size / x_order_ratio
band_y_order    = size / y_order_ratio
band_rcmk_order = size / rcmk_order_ratio

average_band_ratios = [np.mean(no_order_ratio), np.mean(x_order_ratio), np.mean(y_order_ratio), np.mean(rcmk_order_ratio)]
reorder_labels      = ['No reordering', 'X reordering', 'Y reordering', 'RCMK reordering']

# Solvers complexity and band size complexity (constant factors are arbitrary)
complexity_full_solver   = np.power(nodes, 3) / 100000000
complexity_band_solver   = np.power(nodes, 1.65) / 5000000
complexity_no_reordering = size / 1.3
complexity_x_reordering  = nodes ** 0.5 * 13

fig, axes = plt.subplots(2, 1, figsize=(8, 10))

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

plt.tight_layout()
# plt.savefig('../../Rapport/data/solver_plot.pdf', bbox_inches='tight')
plt.show()