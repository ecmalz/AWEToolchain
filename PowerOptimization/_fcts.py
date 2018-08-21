'''
collecting def
Python Version 2.7 / Casadi version 2.4.1
- Author: Elena Malz, Chalmers 2016
'''

import sys
sys.path.append(r"/usr/local/casadi-py27-v3.3.0/")
import casadi as ca
import casadi.tools as ca
import numpy as np
from scipy import interpolate


def Lagrange_poly(x,y):
    "creates an lagrange polynomial through each data point based on x and y data"
    t = ca.SX.sym('t')
    d = x.shape[0]                          # amount of parameters

    poly = 0
    for j in range(d):                  # for all data points ...
        L = y[j]                        # parameter = fct output
        for r in range(d):
            if r != j:
                L *= (t-x[r])/(x[j]-x[r])
        poly+=L
    lfcn = ca.Function('lfcn', [t],[poly])
    return lfcn

def smooth_Lagrange_poly(x,y):
    t    = ca.SX.sym('t')
    d    = len(x)                       # amount of parameters
    tau  = ca.SX.sym('tau',d)              # parameter as minimisation variable
    poly = 0

    for j in range(d):                  # for all data points ...
        L = tau[j]
        for r in range(d):
            if r != j:
                L *= (t-x[r])/(x[j]-x[r])
        poly+=L
    L_fun   = ca.Function('L_fun', [t,tau],[poly])
    ddL,_     = ca.hessian(poly,t)
    ddL_fun = ca.Function('ddL_fun', [t,tau], [ddL])
    # ddL_fun = L_fun.hessian(0)          # second order derivative to
    # [ddL,_,_]  = ddL_fun([t,tau])

    # minimise tau = fct output, incl penalize curvature
    res = 0.1 *  sum([(L_fun(x[k],tau) - y[k])**2 for k in range(d)])[0]
    res += sum([ddL_fun(x[k],tau)[0]**2 * 1e4 for k in range(d)])[0]

    Cost= ca.Function('cost',[tau],[res])
    nlp = {'x': tau, 'f': res}
    solver = ca.nlpsol("solver", "ipopt", nlp)
    sol = solver(**{})
    tau_opt = sol['x']                  # optimal parameter for polynomial
    return  L_fun, tau_opt

def value_Lagrange_poly(x,y,t):
    "creates an lagrange polynomial through each data point based on x and y data and \
    returns directly the value at point t "
    d = x.shape[0]                      # amount of parameters

    poly = 0
    for j in range(d):                  # for all data points ...
        L = y[j]                        # parameter = fct output
        for r in range(d):
            if r != j:
                L *= (t-x[r])/(x[j]-x[r])
        poly+=L
    return poly


def create_LSpoly(d,x,y):
    "creates polynomial by least squares fitting of order d. Builds Vandermonde matrix manually"
    # x = np.append(x,0)              # add (0,0) as data point
    # y = np.append(y,0)

    a = d+1                         # number of parameters including a0
    M = ca.SX.sym('M',2*d+1)           # all the exponents from 1 to 2k
    sizex =  x.shape[0]             # number of data points x
    M[0] = sizex                    # first entry in matrix

    for k in range(1, M.shape[0]): # collect all matrix entries
        M[k] = sum([x[o]**k for o in range(0,sizex)])

    sumM = ca.SX(a,a)                  # create actual matrix
    for j in range(0,a):
            for k in range(0,a):
                sumM[j,k] = M[k+j]

    B = ca.SX(a,1)                     # create B vector (sum of y)
    for k in range(0,a):
        B[k] = sum([y[o]*x[o]**k for o in range(0,sizex)])

    X = ca.solve(sumM, B)    # parameters order: low to high power
    xvar = ca.SX.sym('xvar')
    poly = X[0]
    for k in range(1,X.shape[0]):
        poly += X[k]*xvar**(k)      # create polynomial
    pfun = ca.Function('poly',[xvar],[poly])
    return pfun, X, poly

# x = np.array([2611.55517578,  2350.17895508,  2096.43579102,  1849.80297852,
#         1633.5892334 ,  1469.5637207 ,  1331.72558594,  1196.09936523,
#         1062.48657227,   931.31152344,   802.23907471,   674.7845459 ,
#          548.78527832,   424.55926514,   301.9289856 ,   180.7013855 ,
#           60.22678375, 0.])
# y = np.array([ 700.,  725.,  750.,  775.,  800.,  820.,  835.,  850.,  865.,
#         880.,  895.,  910.,  925.,  940.,  955.,  970.,  985., 988.])    # [hPa] , 1hPa = 100 N / m**2
def rho(h):
    "calculates density of air in heigh h"
    # pressure = interpolate.UnivariateSpline(x, y)
    z = np.array([  4.44105964e-06,  -1.23788651e-01,   9.91302146e+02])
    gasconst = 287.05                 # gas constant
    temp = (8+273)-6.5/1000*h        # T(h) ; 6C as baseline since pressure is from Jan in Got
    return ca.polyval(z,h)*100 / (gasconst*temp)
