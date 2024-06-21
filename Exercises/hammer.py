#
# Computing Hammer quadratures
# Vincent Legat - 2021
# Ecole Polytechnique de Louvain
#

from numpy import *
from numpy.linalg import solve

#
# -1- Hammer : 4 points
#

A = array([[1,3], [1/9,11/25]])
I = array([1/2,1/12])
w = solve(A, I)

print("====== Hammer : 4 points ========= ")
print("              w_1 = %10.7f " % w[0])
print("        w_2 = w_3 = %10.7f " % w[1])                      

#
# -2- Hammer : 4 points
#

a = (6 + sqrt(15)) / 21
b = (6 - sqrt(15)) / 21
A = array([[1, 3, 3],
           [1 / 27, (1 - 2 * b)**3 + 2 * b**3, (1 - 2 * a)**3 + 2 * a**3],
           [1 / 9 , (1 - 2 * b)**2 + 2 * b**2, (1 - 2 * a)**2 + 2 * a**2]])
I = array([1 / 2, 1 / 20, 1 / 12])
w = solve(A, I)

print("====== Hammer : 7 points ========= ")
print("              w_1 = %10.7f " % w[0])
print("  w_2 = w_3 = w_4 = %10.7f " % w[1])                      
print("  w_5 = w_6 = w_7 = %10.7f " % w[2])