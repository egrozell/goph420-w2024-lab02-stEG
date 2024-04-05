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
    Ea = 1.0e-4
    maxit = 30
    eps_a = np.array([])
    it = 0

    for i in range(maxit):
        x1 = x0 - f(x0)/dfdx(x0)
        ea = abs((x1-x0)/x1)
        eps_a = np.append(eps_a, ea)
        it += 1
        if ea < Ea:
            break
        x0 = x1

    return x0, it, eps_a


def root_secant_modified(x0, dx, f):
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

    Ea = 1.0e-7
    maxit = 30
    eps_a = np.array([])
    it = 0

    for i in range(maxit):
        x1 = x0 - f(x0)*dx/(f(x0+dx)-f(x0))
        ea = abs((x1-x0)/x1)
        eps_a = np.append(eps_a, ea)
        it += 1
        if ea < Ea:
            break
        x0 = x1

    return x0, it, eps_a
