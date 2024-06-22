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


def gaussVec(term):
    ids = identifier(term)
    arrays = []
    for i in range(len(ids)):
        if isinstance(eval(ids[i]), np.ndarray):
            arrays.append(ids[i])
    arrayLength = len(eval(arrays[0]))
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    res = []
    for k in range(arrayLength):
        values = []
        for ident in ids:
            if ident in arrays:
                exec("values.append(" + ident + "[k]" + ")")
            else:
                exec("values.append(" + ident + ")")
        derivatives = []
        i = 0
        while i < len(symbols):
            r = sym.diff(termSymbol, symbols[i])
            j = 0
            while j < len(symbols):
                r = r.evalf(subs={symbols[j]: values[j]})
                j += 1
            derivatives.append(r.evalf())
            i += 1
        i = 0
        sigmaArrays = []
        for t in range(len(ids)):
            if isinstance(eval("sigma_" + ids[t]), np.ndarray):
                sigmaArrays.append(ids[t])
        while i < len(derivatives):
            if ids[i] in sigmaArrays:
                exec("derivatives[i] *= sigma_" + ids[i] + "[k]")
            else:
                exec("derivatives[i] *= sigma_" + ids[i])
            i = i + 1
        resj = 0
        for z in derivatives:
            resj += z ** 2
        res.append(sqrt(resj))
    return array(res)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

W1 = matrix("""
1000    340    2   ;
1400    260    2;
1900    350    4;
2500    270    4;
2900    300    5;
3600    290    6;
4000    340    8
""")# Luft; frequenz, delta l/n in mm
#n ist die zahl der gemesenen abstände
#danach verdoplung; lamda=(2*delta l/n)

f1 = toArray(W1[:, 0])
sigma_f1 = 10*ones(len(f1))
delta_L1 = toArray(W1[:, 1])/1000
sigma_delta_L1 = 10/1000*ones(len(delta_L1))
n1 = toArray(W1[:, 2])
sigma_n1 = 0

X1 = n1/(2*delta_L1)
sigma_X1 = gaussVec("n1/(2*delta_L1)")


def linear(x, a, b):
    return a*x + b


errorbar(X1, f1, sigma_f1, sigma_X1,'x', label='Measurement air')
optimizedParameters1, s = opt.curve_fit(linear, X1, f1)
plot(X1, linear(X1, *optimizedParameters1), label="fit")
xlabel('n/2L in 1/m', fontsize=20)
ylabel('Frequency in Hz', fontsize=20)
legend(fontsize=13)
grid()
plt.tight_layout()
savefig('luft')
show()


c1 = optimizedParameters1[0]
sigma_c1 = np.diag(s)[0]

latexTable(UC(f1, 'Hz'), UC(delta_L1, 'mm'), UC(n1))

W2 = matrix("""
1000    180    1;
1500    270    3;
1900    300    4;
2300    250    4;
2800    290    6;
3600    330    8;
4000    290    8
""")# CO2; frequenz, delta l/n in mm
#n ist die zahl der gemesenen abstände
#danach verdoplung; lamda=(2*delta l/n)

T = 24.6 + 273.15


f2 = toArray(W2[:, 0])
sigma_f2 = 10*ones(len(f2))
delta_L2 = toArray(W2[:, 1])/1000
sigma_delta_L2 = 10/1000*ones(len(delta_L2))
n2 = toArray(W2[:, 2])
sigma_n2 = 0

X2 = n2/(2*delta_L2)
sigma_X2 = gaussVec("n2/(2*delta_L2)")


def linear(x, a, b):
    return a*x + b


errorbar(X2, f2, sigma_f2, sigma_X2,'x', label='Measurement CO2')
optimizedParameters2, s = opt.curve_fit(linear, X2, f2)
plot(X2, linear(X2, *optimizedParameters2), label="fit")
xlabel('n/2L in 1/m', fontsize=20)
ylabel('Frequency in Hz', fontsize=20)
legend(fontsize=13)
grid()
plt.tight_layout()
savefig('co2')
show()


c2 = optimizedParameters2[0]
sigma_c2 = np.diag(s)[0]

latexTable(UC(f2, 'Hz'), UC(delta_L2, 'mm'), UC(n2))




# zweiter teil: stab
# einspannung bei L/2
f1 = matrix("""
2150;
6450;
10750;
15050;
19350
""")

L = 113.9/100
sigma_L = 0.5/100

f1 = toArray(f1)
sigma_f1 = 50
n1 = array([1, 3, 5, 7, 9])

errorbar(n1, f1, sigma_f1, None,'x', label='Measurement, fixed at L/2')
optimizedParameters1, s = opt.curve_fit(linear, n1, f1)
plot(n1, linear(n1, *optimizedParameters1), label="fit")
xlabel('Order of the frequency', fontsize=20)
ylabel('Frequencyin Hz', fontsize=20)
legend(fontsize=13)
grid()
plt.tight_layout()
savefig('spann1')
show()


steig = optimizedParameters1[0]
sigma_steig = sigma_f1
c_spann1 = 2*L*steig
sigma_c_spann1 = gauss("2*L*steig")

latexTable(UC(f1, 'Hz'), UC(n1))


#einspannung bei 1/4 3/4
f2 = matrix("""
4300;
8900;
12900;
17500
""")


f2 = toArray(f2)
sigma_f2 = 50
n2 = array([1, 2, 3, 4])

errorbar(n2, f2, sigma_f2, None,'x', label='Measurement, fixed at L/4 and 3L/4')
optimizedParameters2, s = opt.curve_fit(linear, n2, f2)
plot(n2, linear(n2, *optimizedParameters2), label="fit")
xlabel('Order of the frequency', fontsize=20)
ylabel('Frequency in Hz', fontsize=20)
legend(fontsize=13)
grid()
plt.tight_layout()
savefig('spann2')
show()


steig2 = optimizedParameters2[0]
sigma_steig2 = sigma_f2
c_spann2 = L*steig2
sigma_c_spann2 = gauss("2*L*steig2")

latexTable(UC(f2, 'Hz'), UC(n2))

