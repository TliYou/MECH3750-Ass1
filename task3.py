#!/usr/bin/env python

"""@package docstring
MECH3750 Engineering Analysis II
Assignment 1: Buckling analysis

Authors: Merrick Heley & Sophia Impiccini

Main file for Task 3
"""

import math
from prettytable import PrettyTable
import pylab

import integrator as I
from mybisect import solvebisect

limit = 0.99999999

def f1(theta, alpha):
    """
    Evaluate a mathematical function of theta and alpha
    
    Parameters
    ----------

    theta : function
    alpha : float

    Returns
    -------
    
    Value of the function
    """  
        
    return -math.cos(theta)/math.sqrt(math.cos(theta) - math.cos(alpha))

def f2(phi, alpha):
    """
    Evaluate a mathematical function of theta and alpha
    
    Parameters
    ----------

    theta : function
    alpha : float

    Returns
    -------
    
    Value of the function
    """  

    return 1.0/math.sqrt(1.0 - ((math.sin(alpha/2.0))**2 * (math.sin(phi))**2))

def f3(theta, alpha):
    """
    Evaluate a mathematical function of theta and alpha
    
    Parameters
    ----------

    theta : function
    alpha : float

    Returns
    -------
    
    Value of the function
    """
        
    return -math.sin(theta)/math.sqrt(math.cos(theta) - math.cos(alpha)) 

def eqn1(alpha, XgOnL, intg):
    """
    Evaluate the equation for X_g/L given alpha, and subtract the given value
        of XgOnL.
        
    Can be used in optimisers to find a value for alpha.
    
    Parameters
    ----------
    
    alpha : float
        Value of alpha in the function
    XgOnL : float
        Value for X_g/L
    intg : integrator class
        Initialised integrator that can be used to get values
    
    Returns
    -------
    
    Float X_gOnL(alpha) - XgOnL
    """
    
    Term1 = 1.0/(2.0**0.5)
    Term2 = intg.get(f1, alpha, 0, 1e-6, alpha)
    Term3 = intg.get(f2, 0, math.pi/2.0, 1e-6, alpha)
    Term4 = XgOnL
    
    return Term1 * (Term2/Term3) - Term4

def eqn2(alpha, intg):
    """
    Evaluate an equation for Y_b/L
    
    Parameters
    ----------

    alpha : float
        Value of alpha in the function
    intg : integrator class
        Initialised integrator that can be used to get values

    Returns
    -------
    
    Float for Y_b/L
    """  
    
    Term1 = 1.0/(2.0**0.5)
    Term2 = intg.get(f3, alpha, 0, 1e-6, alpha)
    Term3 = intg.get(f2, 0, math.pi/2.0, 1e-6, alpha)
    
    return Term1 * (Term2/Term3)

def eqn3(alpha, intg):
    """
    Evaluate an equation for P/P_e
    
    Parameters
    ----------

    alpha : float
        Value of alpha in the function
    intg : integrator class
        Initialised integrator that can be used to get values

    Returns
    -------
    
    Float for P/P_e
    """  
    
    Term1 = 4.0/(math.pi**2)
    Term2 = intg.get(f2, 0, math.pi/2.0, 1e-6, alpha)**2
    
    return Term1 * Term2

if __name__ == '__main__':
    print("Solve for alpha")
    
    intg = I.Integrator(I.gaussian2, I.gaussian4)
    
    ListA = []
    ListX = []
    ListY = []
    ListP = []
    
    table = PrettyTable(["X_g/L", "alpha", "Y_b/L", "P/P_e"])
    
    for XgOnL in (0.99, 0.95, 0.9, 0.5, 0):
        ListX.append(XgOnL)
        ListA.append(solvebisect(eqn1, 0.01, limit*math.pi, 1e-6, XgOnL, intg))
        ListY.append(eqn2(ListA[-1], intg))
        ListP.append(eqn3(ListA[-1], intg))
        
        table.add_row([XgOnL, math.degrees(ListA[-1]), ListY[-1], ListP[-1]])
    
    print table
    
    pylab.title("Buckling Analysis")
    pylab.xlabel("alpha (radians)"); pylab.ylabel("ratio")
    Plot1, = pylab.plot(ListA, ListX, 'b-')
    Plot2, = pylab.plot(ListA, ListY, 'r-')
    Plot3, = pylab.plot(ListA, ListP, 'g-')
    pylab.legend([Plot1, Plot2, Plot3], ["x_g/L", "y_b/L", "P/Pe"], loc=2)
    #pylab.show()
    pylab.savefig('graph.png')