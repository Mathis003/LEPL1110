import numpy as np
import matplotlib.pyplot as plt

def f(x, A, b, c):
    return float(0.5 * x.T * A * x - b.T * x + c)

def contoursteps(x1, x2, zs, steps=None):
    fig = plt.figure(figsize=(6, 6))
    cp = plt.contour(x1, x2, zs, 10)
    plt.clabel(cp, inline=1, fontsize=10)
    if steps is not None:
        steps = np.matrix(steps)
        plt.plot(steps[:, 0], steps[:, 1], '-o')
    fig.show()


A = np.matrix([[3.0, 2.0], [2.0, 6.0]])
b = np.matrix([[2.0], [-8.0]])  # We will use the convention that a vector is a column vector
c = 0.0
    
x = np.matrix([[-2.0], [2.0]])
steps = [(-2.0, 2.0)]
i = 0
imax = 20
eps = 0.01
r = b - A * x
delta = r.T * r
delta0 = delta

while (i < imax and delta > eps**2 * delta0):
    alpha = 2 * float(delta / (r.T * (A * r)))
    x += alpha * r
    steps.append((x[0, 0], x[1, 0]))  # Store steps for future drawing
    r = b - A * x
    delta = r.T * r
    i += 1

size = 30
zs = np.zeros((size, size))
x1 = list(np.linspace(-6, 6, size))
x2 = list(np.linspace(-6, 6, size))
x1, x2 = np.meshgrid(x1, x2)
for i in range(size):
    for j in range(size):
        x = np.matrix([[x1[i, j]], [x2[i, j]]])
        zs[i, j] = f(x, A, b, c)

contoursteps(x1, x2, zs, steps)

plt.show()