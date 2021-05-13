from scipy.optimize import fsolve
import vector_functions as vt
import Stats_functions as st
import partial_differentiation_functions as pde
import math
import sympy as sym
import numpy as np
from numpy import linalg as la
from scipy import integrate
from sympy import cos, sin, pi
from sympy.vector import ParametricRegion, vector_integrate
from sympy import *
init_printing()

x, y, z, i, j, k, t, u, v, a, b, c, r, theta, phi = sym.symbols('x y z i j k t u v a b c r theta phi')


function_1 = [a*r*cos(phi)*sin(theta),b*r*sin(phi)*sin(theta),c*r*cos(theta)]
function_2 = [x-z,z**2,z]
function_3 = (phi)
order = [x,y,z]
x_bounds = [0,pi]
y_bounds = [0,a]
z_bounds = [0,3]
#p = [1,-3,2]



wrt1 = vt.vector_partial_derivative(function_1,theta)
wrt2 = vt.vector_partial_derivative(function_1,phi)
#pprint(wrt1)
#pprint(wrt2)

cross = vt.cross_product(wrt1,wrt2)

#pprint(cross)

modulus = vt.magnitudealternate(cross)

#mod_squared = modulus**2

#pprint(sim)

pprint(vt.triple_matrix_det(([sin(u),2,0],[7,1,2],[0,5,1])))

X = a*r*sin(theta)*cos(phi)
Y = b*r*sin(theta)*sin(phi)
Z = c*r*cos(theta)
term = [X, Y, Z]
variables = [r, theta, phi]
jacobian = []
for i in range(3):
   new_row = []
   differential = term[i]
   for j in range(3):
       value = sym.diff(differential, variables[j])
       new_row.append(value)
   jacobian.append(new_row)
array = np.array(jacobian)
det = vt.triple_matrix_det(array)
det1 = sym.simplify(det)
pprint(det1)