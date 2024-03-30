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

    zeta_max = np.sqrt(H**2 * (beta1**-2-beta2**-2))

    def initial(frq, _k):
        return ((2*_k+1)/(4*frq)-zeta_max*1e-4)

    freq_list = np.array(np.linspace(1, 100, 20))
    cl_mode = []
    wvl_mode = []
    zeta_mode = []

    for j, freq in enumerate(freq_list):
        zeta = []
        cl = []
        wvl = []
        k = 0
        x0 = initial(freq, k)
        while x0 < zeta_max:

            def g(zeta):
                return (((rho2/rho1)*((np.sqrt(zeta_max**2-zeta**2))/zeta))
                        - (np.tan(2*np.pi*freq*zeta)))

            def dgdx(zeta):
                return ((-rho2*zeta_max**2/(rho1*zeta**2 *
                                            np.sqrt(zeta_max**2-zeta**2)))
                        - (2*np.pi*freq/(np.cos(2*np.pi*freq*zeta)**2)))

            zeta_k, _, _ = root_newton_raphson(x0, g(x0), dgdx(x0))
            zeta.append(zeta_k)
            cl.append(1/np.sqrt(beta1**-2 - (zeta_k/H**2)))
            wvl.append(cl[-1]/freq)

            k += 1
            x0 = initial(freq, k)
        cl_mode.append(cl)
        wvl_mode.append(wvl)
        zeta_mode.append(zeta)

if __name__ == "__main__":
    main()
