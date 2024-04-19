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
                a = plt.pcolormesh(x, y, f, shading="gouraud", linewidth=1, edgecolors="k", *args, vmin=vmin, vmax=vmax, **kwargs)
        return a
    
if __name__ == "__main__":

    nameMesh = "../../Rapport/data/mesh.txt"
    mesh = Mesh(nameMesh)
    nameConstrainsFile = "../../Rapport/data/maxConstrainBridge.txt"
    maxConstrainSol = np.loadtxt(nameConstrainsFile)

    nameSolution = "../../Rapport/data/UV.txt"
    uv = np.loadtxt(nameSolution, skiprows=1, delimiter=",")
    uv_norm = np.linalg.norm(uv, axis=1)

    factor = 2e4

    cb = mesh.plotfield(maxConstrainSol, uv * factor, cmap="cool")
    plt.colorbar(cb)
    mesh.plot(uv * factor, lw=0.2, c="k")
    plt.gca().set_aspect("equal")
    plt.grid(alpha=0.2)
    plt.title('Analysis of Bridge Plastic Deformation:\nVon Mises Equivalent Criterion per Node')
    # plt.savefig("../../Rapport/data/maxConstrainBridge.pdf", bbox_inches='tight')
    plt.show()

# %%