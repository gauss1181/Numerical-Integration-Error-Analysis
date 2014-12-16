# An Analysis of Error Estimates in Numerical Integration
# Megan Chen
# 21-660 Numerical Analysis

# This is all the code I used to calculate numerical approximations
# for each quadrature rule I analyzed in my paper.

from math import *

# Functions/Integrands

def f1(x):
    return 1.0/(5-4*cos(x))

def f2(x):
    return 1.0/((1-x**2)**0.5)

def f3(x):
    return x*sin(x)*cos(x)

def f4(x):
    return x**2-2*x+3

def f5(x):
    return x/((1-x**2)**0.5)

# Trapezoidal Rule

def trap_f1(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [0.5*f1(2*pi*k/n)]
        else:
            terms += [f1(2*pi*k/n)]
    return ((2*pi)/n)*sum(terms)

# Trapezoidal Rule won't work on f2 since it's undefined on the
# integration interval boundaries, x=-1 and x=1.

def trap_f3(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [0.5*f3(2*pi*k/n)]
        else:
            terms += [f3(2*pi*k/n)]
    return ((2*pi)/n)*sum(terms)

def trap_f4(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [0.5*f4(-1+2*k/float(n))]
        else:
            terms += [f4(-1+2*k/float(n))]
    return (2.0/n)*sum(terms)

# Trapezoidal Rule won't work on f5 since it's undefined on the
# integration interval boundaries, x=-1 and x=1.

# Simpson's Rule

def simpson_f1(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [1]
        elif (k % 2 == 1):
            terms += [4*f1(2*pi*k/n)]
        else:
            terms += [2*f1(2*pi*k/n)]
    return ((2*pi)/(3*n))*sum(terms)

# Simpson's Rule won't work on f2 since it's undefined on the
# integration interval boundaries, x=-1 and x=1.

def simpson_f3(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [0]
        elif (k % 2 == 1):
            terms += [4*f3(2*pi*k/n)]
        else:
            terms += [2*f3(2*pi*k/n)]
    return ((2*pi)/(3*n))*sum(terms)

def simpson_f4(n):
    terms = []
    for k in xrange(n+1):
        if (k == 0 or k == n):
            terms += [f4(-1+2*k/float(n))]
        elif (k % 2 == 1):
            terms += [4*f4(-1+2*k/float(n))]
        else:
            terms += [2*f4(-1+2*k/float(n))]
    return (2.0/(3*n))*sum(terms)

# Simpson's Rule won't work on f5 since it's undefined on the
# integration interval boundaries, x=-1 and x=1.

# Gauss-Chebyshev Quadrature (Only works when integrand has
# something divided by sqrt(1-x^2))

def gaussChebyshev_f2(n):
    terms = []
    for k in xrange(n):
        terms += [1]
    return (pi/n)*sum(terms)

def gaussChebyshev_f5(n):
    terms = []
    for k in xrange(n):
        terms += [cos(((2*k+1)*pi)/(2*n))]
    return (pi/n)*sum(terms)

# Gauss-Legendre Quadrature: First Formula (Only works when
# interval of integration is [-1, 1]

def gaussLegendre1_f2():
    return 2*f2(0)

def gaussLegendre1_f4():
    return 2*f4(0)

def gaussLegendre1_f5():
    return 2*f5(0)

# Gauss-Legendre Quadrature: Second Formula (Only works when
# interval of integration is [-1, 1]

def gaussLegendre2_f2():
    return f2(-1.0/3**0.5)+f2(1.0/3**0.5)

def gaussLegendre2_f4():
    return f4(-1.0/3**0.5)+f4(1.0/3**0.5)

def gaussLegendre2_f5():
    return f5(-1.0/3**0.5)+f5(1.0/3**0.5)

# Rectangular Rule for 2*pi-Periodic Functions
# (Only works for 2*pi-Periodic integrands)

def rectangular_f1(n):
    terms = []
    for k in xrange(1, n+1):
        terms += [f1(2*pi*k/n)]
    return (2*pi/n)*sum(terms)

def rectangular_f3(n):
    terms = []
    for k in xrange(1, n+1):
        terms += [f3(2*pi*k/n)]
    return (2*pi/n)*sum(terms)
