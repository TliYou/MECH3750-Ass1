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
    Adaptive integrator class that is initialised with two integration rules.
    """
    
    def __init__(self, integrator1, integrator2):
        """ Constructor. 
        
        Initialises with two integrator functions of the form:
            integrate(f, a, b, *args)
        
            Parameters
            ----------
        
            func : function
                A Python function or method to integrate.
            a : float
                Lower limit of integration.
            b : float
                Upper limit of integration.
            args : tuple, optional
                extra arguments to pass to func
        
            Returns
            -------
            
            value of the integral
        """
        
        self.intg1 = integrator1
        self.intg2 = integrator2
        
    def get(self, f, a, b, tol, *args):
        """ 
        Get the integral of a function between a and b to accuracy of tol.
        Returns the numeric integral.
        
        Parameters
        ----------
    
        func : function
            A Python function or method to integrate.
        a : float
            Lower limit of integration.
        b : float
            Upper limit of integration.
        tol : float
            Acceptable tolerance for the difference between integration funcs
        args : tuple, optional
            extra arguments to pass to func
    
        Returns
        -------
        
        value of the integral
        """
        
        I1 = self.intg1(f, a, b, *args)
        I2 = self.intg2(f, a, b, *args)
        if abs(I2 - I1) > tol:
            mid = 0.5 * (a + b)
            I = self.get(f, a, mid, tol, *args) \
                    + self.get(f, mid, b, tol, *args)
        else:
            I = I2
        
        return I

def _gaussian(f, a, b, w, x, *args):
    """
    Calculate an integral of f between a and b
        using the given weights and values of a gaussian quadrature function
        
    Parameters
    ----------

    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    w : list
        List of weightings for quadrature
    x : list
        List of corresponding points for quadrature
    args : tuple, optional
        extra arguments to pass to func

    Returns
    -------
    
    value of the integral    
    """
    
    fv = [(b-a)/2.0*x[i] + (a+b)/2.0 for i in xrange(len(x))]
    return ((b-a)/2.0)*sum(w[i]*f(fv[i], *args) for i in xrange(len(w)))

def gaussian2(f, a, b, *args):
    """
    Calculate an integral of f between a and b
        using 2-pt gaussian quadrature
        
    Parameters
    ----------

    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    args : tuple, optional
        extra arguments to pass to func

    Returns
    -------
    
    value of the integral
    """
    
    w = [1, 1]
    x = [-math.sqrt(3)/3.0, math.sqrt(3)/3.0]
    return _gaussian(f, a, b, w, x, *args)

def gaussian4(f, a, b, *args):
    """
    Calculate an integral of f between a and b
        using 4-pt gaussian quadrature
        
    Parameters
    ----------

    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    args : tuple, optional
        extra arguments to pass to func

    Returns
    -------
    
    value of the integral
    """
    
    w = [(18-math.sqrt(30))/36.0, (18+math.sqrt(30))/36.0, 
         (18+math.sqrt(30))/36.0, (18-math.sqrt(30))/36.0]
    x = [-math.sqrt((3+2*math.sqrt(6.0/5.0))/7.0), 
         -math.sqrt((3-2*math.sqrt(6.0/5.0))/7.0), 
         math.sqrt((3-2*math.sqrt(6.0/5.0))/7.0), 
         math.sqrt((3+2*math.sqrt(6.0/5.0))/7.0)]
    return _gaussian(f, a, b, w, x, *args)
        
if __name__ == '__main__':
    print("Test Adaptive Integrator")
    
    f = [lambda x: 5*x + 3,
         lambda x: x**3 + 3**x + 5*x**2,
         lambda x: 5*x**5,
         lambda x: 5*x**5 - x**4 + 3**x + math.sqrt(abs(x)) + x*math.sin(x**2)]
    
    integrator = Integrator(gaussian2, gaussian4)
    
    for i in f:
        print "Scipy Integral: {}".format(scipy.integrate.quad(i, -5, 5)[0])
        print "My Integral: {}".format(integrator.get(i, -5, 5, 1e-6))
        print "----"