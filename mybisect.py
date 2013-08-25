#!/usr/bin/env python

"""@package docstring
MECH3750 Engineering Analysis II
Assignment 1: Buckling analysis

Authors: Merrick Heley & Sophia Impiccini

Includes bisection solver
"""

import scipy.optimize as opt
import math

def solvebisect(f, a, b, tol, *args):
    """
    Use bisection method to find roots between a and b of function f, 
        with an accuracy of tol.
        
    Parameters
    ----------

    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    tol : float
        Acceptable tolerance for the difference between points
    args : tuple, optional
        extra arguments to pass to func

    Returns
    -------
    
    Point of the root    
        
    """

    fa = f(a, *args)
    fb = f(b, *args)
    c = 0.5*(a+b)

    if abs(fa) < tol:   # Check if root is at extrema
        return a

    if abs(fb) < tol:   # Check if root is at extrema
        return b

    if fa*fb >= 0.00:   # Check if root within range
        return None

    while True: #abs(b-a)/2.0 > tol:   #Run while outside tolerance
        c = 0.5*(a+b)
        fa, fc = f(a, *args), f(c, *args)

        if abs(fc) < tol:	# Quit when within tolerance
            return c

        if fa*fc < 0.0:	# Recalculate a or b
            b = c
        else:
            a = c
    return c

if __name__ == '__main__':
    print("Test Bisection Method")
    f = [lambda x: 5*x + 3,
         lambda x: x**3 + 3**x + 5*x**2,
         lambda x: 5*x**5,
         lambda x: 5*x**5 - x**4 + 3**x + math.sqrt(abs(x)) + x*math.sin(x**2)]

    for i in f:
        print "Scipy Solution: {}".format(opt.bisect(i, -10, 10, maxiter=200))
        print "My Solution: {}".format(solvebisect(i, -10, 10, 1e-12))
        print "----"
