import numpy as np
import matplotlib.plot as plt
from goph420_lab02.root_finder import (root_secant_modified,
                                       root_newton_raphson)

def main():

    rho1 = 1800.0  # kg/m**3
    rho2 = 2500.0  # kg/m**3
    beta1 = 1900.0  # m/s
    beta2 = 3200.0  # m/s
    H = 4.0  # km
    eps_t = 1.0e-8

    zeta_max = np.sqrt(H**2 * (beta1**-2-beta2**-2))

    def g(zeta, freq):
        return (((rho2/rho1)*((np.sqrt(zeta_max**2-zeta**2))/zeta))
                - (np.tan(2*np.pi*freq*zeta)))

    def dgdx(zeta, freq):
        return ((-rho2*zeta_max**2/(rho1*zeta**2*np.sqrt(zeta_max**2-zeta**2)))
                - (2*np.pi*freq/(np.cos(2*np.pi*freq*zeta)**2)))

    freq_list = []
    freq_mode = []
    cl_mode = []
    lambda_mode = []



if __name__ == "__main__":
    main()
