#
# Computing P2-C0 shape functions for triangles
#
# Vincent Legat - 2021
# Ecole Polytechnique de Louvain
#

from numpy import *
from numpy.linalg import solve
import matplotlib 
import matplotlib.pyplot as plt 
matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet

# =========================================================================

def shapeComputeCoefficient(nodes):
  n = len(nodes)
  B = diag(ones(n))
  A = zeros([n,n])
  for i in range(0,n):
    xi  = nodes[i,0]
    eta = nodes[i,1]
    A[i,:] = [1, xi, eta, xi*xi, xi*eta, eta*eta, xi**3, xi*xi*eta, xi*eta*eta, eta**3 ]
  return transpose(solve(A,B))
  
# =========================================================================

def shapeCompute(phi,xi,eta):
  return (phi[0] + phi[1]*xi + phi[2]*eta 
                 + phi[3]*xi**2 + phi[4]*xi*eta + phi[5]*eta**2 
                 + phi[6]*xi**3 + phi[7]*eta*xi**2 + phi[8]*xi*eta**2 
                 + phi[9]*eta**3);

# =========================================================================

#
# -1- Computing coefficients
#

nodes = array([[0,0],[1,0],[0,1],
              [1/3,0],[2/3,0],[2/3,1/3],
              [1/3,2/3],[0,2/3],[0,1/3],[1/3,1/3]])
n = len(nodes)              
coeff = shapeComputeCoefficient(nodes)
print(" ==== Coefficients ===== ")
print(coeff)

#
# -2- Checks if the sum is equat to one 
#     for a random point (0.2;0.4)
#

sumPhi = 0;
for i in range(0,n):
  X = 0.2; Y = 0.4
  sumPhi += shapeCompute(coeff[i,:],X,Y)
print(" ==== Check of the sum of shape fuctions : %14.7e ===== " % sumPhi)

#
# -4 Nice plot
#

plt.figure("P2-C0 Shapes functions")
for i in range(0,n): 
  Xg,Yg = meshgrid(linspace(0,1,100),linspace(0,1,100))
  Xt = [0,0,1,0];   Yt = [0,1,0,0];  
  U = shapeCompute(coeff[i,:],Xg,Yg)
  U[Xg+Yg>1] = nan
  
  Xg = Xg + 3.5*nodes[i,0];
  Yg = Yg + 3.5*nodes[i,1];
  Xt = Xt + 3.5*nodes[i,0];
  Yt = Yt + 3.5*nodes[i,1];
  plt.contourf(Xg,Yg,U,10,cmap=myColorMap,vmin=-0.3,vmax=1)
  plt.contour(Xg,Yg,U,10,colors='k',linewidths=1)
  plt.plot(Xt,Yt,'-k');
  
plt.axis("equal"); plt.axis("off")
plt.show()