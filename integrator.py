#!/usr/bin/env python

"""@package docstring
MECH3750 Engineering Analysis II
Assignment 1: Buckling analysis

Authors: Merrick Heley & Sophia Impiccini

Includes an adaptive integrator based on Gaussian quadrature.
"""

import scipy.integrate
import math

class Integrator:
    """
    Adaptive integrator
    """
    
    def __init__(self, integrator1, integrator2):
        """ Constructor. 
        Takes 2 integrator methods of the form intg1(f, a, b)        
        """
        self.intg1 = integrator1
        self.intg2 = integrator2
        
    def get(self, f, a, b, tol=1e-6):
        """ 
        Get the integral of a function between a and b to accuracy of tol.
        Returns the numeric integral.
        """
        I1 = self.intg1(f, a, b)
        I2 = self.intg2(f, a, b)
        if abs(I2 - I1) > tol:
            mid = 0.5 * (a + b)
            I = self.get(f, a, mid, tol) + self.get(f, mid, b, tol)
        else:
            I = I2
        
        return I

def gaussian(f, a, b, w, x):
    """
    Calculate an integral of f between a and b
        using the given weights and values of a gaussian quadrature function
    """
    fv = [(b-a)/2.0*x[i] + (a+b)/2.0 for i in xrange(len(x))]
    return ((b-a)/2.0)*sum(w[i]*f(fv[i]) for i in xrange(len(w)))

def gaussian2(f, a, b):
    """
    Calculate an integral of f between a and b
        using 2-pt gaussian quadrature
    """
    w = [1, 1]
    x = [-0.57735027, 0.57735027]
    return gaussian(f, a, b, w, x)

def gaussian4(f, a, b):
    """
    Calculate an integral of f between a and b
        using 4-pt gaussian quadrature
    """
    w = [0.3478548, 0.6521451, 0.6521451, 0.3478548]
    x = [-0.86113631, -0.33998104, 0.33998104, 0.86113631]
    return gaussian(f, a, b, w, x)
        
if __name__ == '__main__':
    print("Test Adaptive Integrator")
    
    f = [lambda x: 5*x + 3,
         lambda x: x**3 + 3**x + 5*x**2,
         lambda x: 5*x**5,
         lambda x: 5*x**5 - x**4 + 3**x + math.sqrt(abs(x)) + x*math.sin(x**2)]
    
    integrator = Integrator(gaussian2, gaussian4)
    
    for i in f:
        print "Scipy Integral: {}".format(scipy.integrate.quad(i, -5, 5)[0])
        print "My Integral: {}".format(integrator.get(i, -5, 5))
        print "----"