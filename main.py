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
    
def trap(f, a, b):
    """
    Calculate an integral of f between a and b
        using trapezoidal rule (2-pt Newton-Cotes)
    """
    return ((b-a)/2.0)*(f(a) + f(b))

def simpsons4(f, a, b):
    """
     Calculate an integral of f between a and b
         using simpsons 3/8th's method (4-pt Newton-Cotes)
     """
    return ((b-a)/8.0)*(f(a) + 3.0*f((2*a+b)/3.0) + 3.0*f((a+2*b)/3.0) + f(b))
        
if __name__ == '__main__':
    print("Test Adaptive Integrator")
    
    f = [lambda x: 5*x + 3,
         lambda x: x**3 + 3**x + 5*x**2,
         lambda x: 5*x**5,
         lambda x: 5*x**5 - x**4 + 3**x + math.sqrt(abs(x)) + x*math.sin(x**2)]
    
    integrator = Integrator(trap, simpsons4)
    
    for i in f:
        print "Scipy Integral: {}".format(scipy.integrate.quad(i, -5, 5)[0])
        print "My Integral: {}".format(integrator.get(i, -5, 5))
        print "----"