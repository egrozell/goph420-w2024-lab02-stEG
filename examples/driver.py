import numpy as np
import matplotlib.plot as plt
from goph420_lab02.root_finder import root_newton_raphson

def main():

    rho1 = 1800.0  # kg/m**3
    rho2 = 2500.0  # kg/m**3
    beta1 = 1900.0  # m/s
    beta2 = 3200.0  # m/s
    H = 4.0  # km

    zeta_max = np.sqrt(H**2 * (beta1**-2-beta2**-2))  # maximum value for zeta

    # initial root guess function
    def initial(frq, _k):
        return ((2*_k+1)/(4*frq)-zeta_max*1e-4)

    # list of frequencies in Hz to calculate initial guess at
    freq_list = np.array(np.linspace(1, 100, 20))
    cl = []  # velocity list which will hold lists for the modes found at freq
    wvl = []  # wavelength list
    zeta = []  # zeta list

    # Iterating over the frequency list
    for j, freq in enumerate(freq_list):
        zeta_mode = []
        cl_mode = []
        wvl_mode = []
        k = 0  # mode number
        x0 = initial(freq, k)  # guess of root for first mode
        while x0 < zeta_max:  # iterating over the modes to get to zeta max

            def g(zeta):
                return (((rho2/rho1)*((np.sqrt(zeta_max**2-zeta**2))/zeta))
                        - (np.tan(2*np.pi*freq*zeta)))

            def dgdx(zeta):
                return ((-rho2*zeta_max**2/(rho1*zeta**2 *
                                            np.sqrt(zeta_max**2-zeta**2)))
                        - (2*np.pi*freq/(np.cos(2*np.pi*freq*zeta)**2)))

            # Using the newton raphson root finding method
            zeta_k, _, _ = root_newton_raphson(x0, g(x0), dgdx(x0))

            # putting data of root for zeta's value
            # and the values that can then be derived from this value
            # into their respective lists
            zeta_mode.append(zeta_k)
            cl_mode.append(1/np.sqrt(beta1**-2 - (zeta_k/H**2)))
            wvl_mode.append(cl_mode[-1]/freq)

            k += 1
            x0 = initial(freq, k)
        # putting the lists for data at the specific frequency in their list
        cl.append(cl_mode)
        wvl.append(wvl_mode)
        zeta.append(zeta_mode)

if __name__ == "__main__":
    main()
