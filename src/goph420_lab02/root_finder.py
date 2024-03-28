import numpy as np 
def root_newton_raphson(x0, f, dfdx):
    """
    inputs
    ______
    x0: float the initial guess of root 
    f: callable function
    dfdx: callable function which is the first derivative of f
    outputs
    _______
    double: float
    integer: int 
    x_array: numpy.ndarray

    raises
    ______ 
    ValueError: if x0 is not float like 
    TypeError: if f is not callable
    TypeError: if dfdx is not callable
    """
    pass
def root_secant_modified(x0, dx ,f):
    """
    inputs
    ______
    x0: float the initial guess of root 
    dx: is a step size for derivative estimation
    f: callable function
    outputs
    _______
    double: float
    integer: int 
    x_array: numpy.ndarray

    raises
    ______ 
    ValueError: if x0 is not float like 
    ValueError: if dx is not float like 
    TypeError: if f is not callable
    """
    pass
