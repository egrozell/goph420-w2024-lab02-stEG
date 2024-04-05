import numpy as np
import matplotlib.pyplot as plt
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
        print(f"freq {frq} mode {_k}")
        return (((2*_k+1)/(4*frq))-1e-2)

    # list of frequencies in Hz to calculate initial guess at
    freq_list = np.array(np.linspace(1, 100, 10))
    print(freq_list)
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
        # print(f"{x0} < {zeta_max}")
        while x0 < zeta_max:  # iterating over the modes to get to zeta max
        # print(x0<zeta_max)
        # for k in freq_list:

            def g(_zeta):
                return ((rho2/rho1)*np.sqrt(zeta_max**2 - _zeta**2)/_zeta - (np.tan(2*np.pi*freq*_zeta)))

            tempg = g(x0)
            print(tempg)

            def dgdx(_zeta):
                return ((-rho2*zeta_max**2/(rho1*_zeta**2 *
                                            np.sqrt(zeta_max**2-_zeta**2))
                        - (2*np.pi*freq/(np.cos(2*np.pi*freq*_zeta)**2))))

            tempgprime =dgdx(x0)
            print(tempgprime)

            # Using the newton raphson root finding method
            zeta_k, itr, _ = root_newton_raphson(x0, g, dgdx)
            if itr > 20:
                break
            # putting data of root for zeta's value
            # and the values that can then be derived from this value
            # into their respective lists
            zeta_mode.append(zeta_k)
            # print(zeta_k)
            cl_mode.append(1/np.sqrt(beta1**-2 - (zeta_k/H**2)))
            wvl_mode.append(cl_mode[-1]/freq)

            k += 1
            x0 = initial(freq, k)
        # putting the lists for data at the specific frequency in their list
        cl.append(cl_mode)
        wvl.append(wvl_mode)
        zeta.append(zeta_mode)

    print("freq_list")
    print(freq_list)
    print()

    print("cl")
    print(cl)
    print()

    # plt.figure(figsize = (6,8))
    #
    # plt.subplot(2,1,1)
    # for f,c in zip(freq_list,cl):
    #     plt.plot(f,c)
    # plt.xlabel('f [Hz]')
    # plt.ylabel('c_L [m/s]')
    # plt.title('figure 1: c_L as a function of f', fontsize = 7)
    # plt.legend([f'mode {k}' for k,_ in enumerate (f)])
    #
    # plt.subplot(2,1,2)
    # for f,lam in zip(freq_list,lambda_mode):
    #     plt.plot(freq_list,lam)
    # plt.xlabel('f [Hz]')
    # plt.ylabel('lambda [m]')
    # plt.title('figure 2: lambda_L as a function of f', fontsize = 7)
    # plt.legend([f'mode {k}' for k,_ in enumerate (freq_list)])  
    #
    # plt.savefig("figures/mode")

if __name__ == "__main__":
    main()
