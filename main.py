import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp, arcsin, arccos, arctan2, sinh, cosh
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)







def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

W1 = matrix("""
1000    340    2   
1400    260    2
1900    350    4
2500    270    4
2900    300    5
3600    290    6
4000    340    8
""")# Luft; frequenz, delta l/n in mm
#n ist die zahl der gemesenen abstände
#danach verdoplung; lamda=(2*delta l/n)

W2 = matrix("""
1000    180    1
1500    270    3
1900    300    4
2300    250    4
2800    290    6
3600    330    8
4000    290    8
""")# CO2; frequenz, delta l/n in mm
#n ist die zahl der gemesenen abstände
#danach verdoplung; lamda=(2*delta l/n)

T = 24.6 + 273.15





# zweiter teil: stab
# einspannung bei L/2
f1 = matrix("""
2150
6450

10750
15050
19350



""")

L = 113.9/100

#einspannung bei 1/4 3/4
f2 = matrix("""
4300
8900
12900
17500
""")


