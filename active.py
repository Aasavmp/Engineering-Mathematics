from scipy.optimize import fsolve
import Stats_functions as st
from vector_functions import *
import vector_functions as vt
from Stats_functions import *
import partial_differentiation_functions as pdf
from partial_differentiation_functions import *
import math
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy import integrate
from sympy import cos, sin, pi
from sympy.vector import ParametricRegion, vector_integrate
from sympy import *
init_printing()

x, y, z, i, j, k, t, u, v, a, b, c, r, s, w, theta, phi, omega = sym.symbols('x y z i j k t u v a b c r s w theta phi omega')

# pprint(vt.flux_integral([(x**3)*y, -(x**2)*z**2, -(x**2)*y*z], [u*cos(v), u*sin(v), sqrt(u**2-1)], [1, 5], [0, 2*pi]))

# pprint(pdf.inverse_laplace_transform((-2*s-2)/(s**2+4*s+8)))

#pprint(pdf.initial_value_problems_laplace([1, 8, 25], cos(3*t), 2, 1))

#pprint(vt.gauss([(x*y**2), (y*cos(z) + y*z**2), (x**2*z-sin(z))], [0, 1], [0, 2*pi], [0, pi], [r*cos(u)*sin(v), r*sin(u)*sin(v), r*cos(v)]))

#pprint(pdf.continuous_fourier_series(t, 2))

function = lambda x: (2*(-1)**(x+1))/((x*pi)) if x != 0 else 0

#function = lambda x: sin(x) if sin(x) != 0 else 0

#pprint(Stats_functions.t_statistic(5.006, 5.105, sqrt(0.124), 50))


#pprint(vt.flux_integral([x*y, x*y, z], [u*cos(v), u*sin(v), u], [0, 1], [0, 2*pi], direction = 'negative'))

#pprint(pdf.initial_value_problems_laplace([1, 2, 10], 10*exp(-2*t), 2, 3))

#pprint(flux_integral([x*exp(y), y*exp(x), x**2], [u*cos(v), u*sin(v), 0], [0, 1], [0, 2*pi]))

#pprint(triple_integral(r, [r, theta, z], [0, 1], [0, 2*pi], [-1, 1]))

#pprint(pdf.half_range_sine_series_continuous(t/pi, pi, t))

# pprint(Stats_functions.goodness_of_fit([873, 77, 50], [821, 162, 17], 0.05, 1))

v = [(-y)/(x**2+y**2), (x)/(x**2+y**2), 0]

pprint(divergent(v))
pprint((curl(v)))