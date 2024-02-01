#
# Steepest descent and conjugate gradients 
# Vincent Legat - 2021
# Ecole Polytechnique de Louvain
#
# largement inspiré du programme de JFR :-)
#

from numpy import *

#
# -1- The problems :-)
#

A = array([[3.0,2.0],[2.0,6.0]])
b = array([2.0,-8.0])
xinit = array([-4,3.0])
Xmin = -4.5; Xmax = 3.5; Ymin = -2.5; Ymax = 3.5

 
# A = array([[2.0,-1.0],[-1.0,2.0]])
# b = array([1.0,1.0])
# xinit = array([1.5,2.0])
# Xmin = 0; Xmax = 2; Ymin = 0.7; Ymax = 2.2

eps = 0.0001;

#
# -2- Steepest descent
#

x = xinit.copy()
steepest = [(x[0],x[1])]; i = 0; delta = 1; imax = 20

while i < imax and delta > eps :
  r = b - A @ x
  delta = r @ r
  alpha = delta / (r @ A @ r)
  x = x + alpha * r
  steepest.append((x[0],x[1])) 
  i += 1
  print(" = Iteration %4d : %14.7e" % (i,delta))
  
steepest = array(steepest)

#
# -3- Conjugate gradients
#

x = xinit.copy()
conjugate = [(x[0],x[1])]; i = 0; delta = 1; imax = 2

r = b - A @ x
d = r.copy()
while i < imax and delta > eps :
  delta = r @ r
  s  = A @ d
  alpha = delta / (r @ s)
  x = x + alpha * d
  r = r - alpha * s
  beta = (r @ r) / delta
  d = r + beta * d
  conjugate.append((x[0],x[1])) 
  i += 1
  print(" = Iteration %4d : %14.7e" % (i,delta))

conjugate = array(conjugate)

#
# -4- Nice plots
#

import matplotlib 
import matplotlib.pyplot as plt 
matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet

n = 30;

X,Y = meshgrid(linspace(Xmin,Xmax,n),linspace(Ymin,Ymax,n))
J = zeros((n,n))
for i in range(n):
  for j in range(n):
    x = array([X[i,j],Y[i,j]])
    J[i,j] = (x @ A @ x / 2 - b @ x)
    

plt.figure("Conjugate gradients - steepest descent")
plt.contourf(X,Y,J,12,cmap=myColorMap)
plt.axis('equal'); plt.axis('off')
plt.plot(steepest[:,0], steepest[:,1], '-ow')
plt.plot(conjugate[:,0], conjugate[:,1], '-or')
plt.show()