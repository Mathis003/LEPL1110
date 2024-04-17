import numpy as np
import matplotlib.pyplot as plt

mesh_size, nodes, full_solver, band_solver_no_order, band_solver_x_order, band_solver_y_order, band_solver_rmck_order, no_order_ratio, x_order_ratio, y_order_ratio, rcmk_order_ratio = np.loadtxt('../../Processing/data/solver_type.csv', delimiter=',', unpack=True, skiprows=1)

average_band_ratios = [np.mean(no_order_ratio), np.mean(x_order_ratio), np.mean(y_order_ratio), np.mean(rcmk_order_ratio)]

reorder_labels = ['No reordering', 'X reordering', 'Y reordering', 'RMCK reordering']

plt.figure("Performance of Different Solvers with Different Reordering")
plt.title("Performance of Different Solvers with Different Reordering")
plt.plot(nodes, full_solver, '--o', label='Full Solver')
plt.plot(nodes, band_solver_no_order, '--o', label='Band Solver without reordering')
plt.plot(nodes, band_solver_x_order, '--o', label='Band Solver with X reordering')
plt.plot(nodes, band_solver_y_order, '--o', label='Band Solver with Y reordering')
plt.plot(nodes, band_solver_rmck_order, '--o', label='Band Solver with RMCK reordering')
plt.xlabel('Number of Nodes [/]')
plt.ylabel('Time [s]')
plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper left', fontsize='small')


# Adding average band ratios to the plot as a box in the bottom right corner
bbox_props = dict(boxstyle="round", facecolor="white", alpha=0.5)
plt.text(0.95, 0.05, '\n'.join([f'{label}: {ratio:.2f}' for label, ratio in zip(reorder_labels, average_band_ratios)]),
         transform=plt.gca().transAxes, verticalalignment='bottom', horizontalalignment='right', bbox=bbox_props)
plt.text(0.95, 0.24, 'Average Band Ratios :', transform=plt.gca().transAxes, verticalalignment='bottom', horizontalalignment='right')

# plt.show()
plt.savefig('../../Processing/data/solver_plot.png')
plt.savefig('../../Processing/data/solver_plot.svg')
plt.savefig('../../Processing/data/solver_plot.pdf')