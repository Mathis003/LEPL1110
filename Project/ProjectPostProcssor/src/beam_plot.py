#%%
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt

class Mesh:
    nlocal : int
    nnodes : int
    nelem : int
    nodes : NDArray[np.float64]
    elem : NDArray[np.int32]

    def __init__(self, fname):
        self.fname = fname

        with open(fname, "r") as f:

            ######## Nodes ########
            self.nnodes = int(f.readline().split(" ")[3])
            self.nodes = np.zeros((self.nnodes, 2))
            for i in range(self.nnodes):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                part_coord = parts[1].split()
                self.nodes[i] = [float(val.strip()) for val in part_coord]

            line = f.readline()
            while("triangles" not in line and "quads" not in line):
                line = f.readline()

            ######## Elements ########
            self.nelem = int(line.split(" ")[3])

            if "triangles" in line: self.nlocal = 3
            else:                   self.nlocal = 4

            self.elem = np.zeros((self.nelem, self.nlocal), dtype=np.int32)
            for i in range(self.nelem):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.elem[i] = [int(val.strip()) for val in parts[1].split()][:self.nlocal]

    def __str__(self):
        return f"Mesh: {self.fname}\n├─nnodes: {self.nnodes}\n├─nelem: {self.nelem}\n├─local: {self.nlocal}\n├─nodes: array(shape=({self.nodes.shape[0]},{self.nodes.shape[1]}), dtype=np.float64)\n└─elem: array(shape=({self.elem.shape[0]},{self.elem.shape[1]}), dtype=np.int32)\n"

    def __repr__(self) -> str:
        return self.__str__()

    def plot(self, displace=None, *args, **kwargs):
        coord = self.nodes if displace is None else np.array(displace) + self.nodes
        if self.nlocal == 3:
            a = plt.triplot(coord[:,0], coord[:,1], self.elem, *args, **kwargs)
        else:
            y = coord[self.elem[:, [0, 1, 2, 3, 0, 0]], 1]
            y[:,5] = np.nan
            y = y.ravel()
            x = coord[self.elem[:, [0, 1, 2, 3, 0, 0]], 0]
            x[:,5] = np.nan
            x = x.ravel()
            a = plt.plot(x, y, *args, **kwargs)

        # Paramètres poutre
        L = 8.0
        rho = 7850
        w = (1e7 + rho * 9.81 * L) / L
        E = 2.11e11
        I = 1**4 / 12

        def y(x):
            return (- w * x**2 / (24 * E * I)) * (x**2 - 4 * L * x + 6 * L**2)

        factor_def = 1e2
        x_values = np.linspace(0, L, 1000)
        y_values = y(x_values) * factor_def +  1/2

        plt.plot(x_values, y_values, "r", lw=2, label="Analytical Solution")
        plt.legend()
        return a

    def plotfield(self, field, displace=None, *args, **kwargs):
        coord = self.nodes if displace is None else np.array(displace) + self.nodes
        field = np.array(field)
        if len(field) != self.nnodes and field.size != self.nnodes:
            raise ValueError("Field must be a 1D array with length equal to number of nodes")
        field = field.reshape(-1)
        if self.nlocal == 3:
            a = plt.tripcolor(coord[:,0], coord[:,1], self.elem, field, shading="gouraud", *args, **kwargs)
        else:
            vmin = kwargs.pop("vmin", field.min())
            vmax = kwargs.pop("vmax", field.max())
            for i in range(self.nelem):
                x = coord[self.elem[i, [[0, 1], [3,2]]], 0]
                y = coord[self.elem[i, [[0, 1], [3,2]]], 1]
                f = field[self.elem[i, [[0, 1], [3,2]]]]
                a = plt.pcolormesh(x, y, f, shading="gouraud", linewidth=1, edgecolors="k",  *args, vmin=vmin, vmax=vmax, **kwargs)
        return a
    
if __name__ == "__main__":

    nameMesh = "../../Rapport/data/mesh_beam.txt"
    mesh = Mesh(nameMesh)

    print(mesh)
    nameSolution = "../../Rapport/data/UV_beam.txt"

    uv = np.loadtxt(nameSolution, skiprows=1, delimiter=",")
    uv_norm = np.linalg.norm(uv, axis=1)
    factor_def = 1e2

    cb = mesh.plotfield(uv_norm, uv * factor_def, cmap="turbo")
    plt.colorbar(cb)
    mesh.plot(uv * factor_def, lw=0.2, c="k")
    plt.gca().set_aspect("equal")
    plt.grid(alpha=0.2)
    plt.title('Elastic Deformation of Beam under\nForce Density on Both Decks')
    # plt.savefig("../../Rapport/data/elasticityBeam.pdf", bbox_inches='tight')
    plt.show()
# %%
