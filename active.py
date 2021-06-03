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

function = lambda x: (1/(x**2*pi**2))*cos(x*pi/2) if x!=0 else 0

print(magnitude_spectrum(function))
