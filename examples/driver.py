import numpy as np
import matplotlib.pyplot as plt
from goph420_lab02.root_finder import root_newton_raphson

def main():

    rho1 = 1800.0   # kg/m**3
    rho2 = 2500.0   # kg/m**3
    beta1 = 1900.0  # m/s
    beta2 = 3200.0  # m/s
    H = 4.0e03      # m

    zeta_max = np.sqrt(H**2 * (beta1**-2-beta2**-2))  # maximum value for zeta

    # initial root guess function
    def initial(frq, _k):
        # print(f"freq {frq} mode {_k}")
        return (((2*_k+1)/(4*frq))-1e-6)

    # list of frequencies in Hz to calculate initial guess at
    freq_list = np.array(np.linspace(1, 100, 50))
    # print(freq_list)
    cl = []  # velocity list which will hold lists for the modes found at freq
    wvl = []  # wavelength list
    zeta = []  # zeta list
    freq_indx = []

    # Iterating over the frequency list
    for j, freq in enumerate(freq_list):
        zeta_mode = []
        cl_mode = []
        wvl_mode = []
        freq_mode = []
        k = 0  # mode number
        x0 = initial(freq, k)  # guess of root for first mode

        def g(_zeta):
            return ((rho2/rho1)*np.sqrt(zeta_max**2 - _zeta**2)/_zeta
                    - (np.tan(2*np.pi*freq*_zeta)))

        def dgdx(_zeta):
            return (-((2 * np.pi * freq) / np.cos(2 * freq * np.pi * _zeta)**2
                    + (rho2 / rho1) * np.sqrt((H**2 * (beta1**-2-beta2**-2))
                    - _zeta ** 2) / _zeta ** 2
                    + (rho2 / rho1) / (np.sqrt((H**2 * (beta1**-2-beta2**-2))
                                               - _zeta ** 2))))

        if x0 > zeta_max:
            x0 = zeta_max - 1e-5

        # print(f"{x0} < {zeta_max}")
        # print("inside sqrt")
        # print((H**2 * (beta1**-2-beta2**-2)) - (x0**2-1e-02))

        while x0 < zeta_max:  # iterating over the modes to get to zeta max

            # Using the newton raphson root finding method
            zeta_k, _, _ = root_newton_raphson(x0, g, dgdx)

            # putting data of root for zeta's value
            # and the values that can then be derived from this value
            # into their respective lists
            zeta_mode.append(zeta_k)
            # print(zeta_k)
            cl_mode.append(1/np.sqrt(beta1**-2 - (zeta_k/H**2)))
            wvl_mode.append(cl_mode[-1]/freq)
            freq_mode.append(freq)

            k += 1
            x0 = initial(freq, k)
        # putting the lists for data at the specific frequency in their list
        cl.append(cl_mode)
        wvl.append(wvl_mode)
        zeta.append(zeta_mode)
        freq_indx.append(freq_mode)

    # print(zeta)
    plt.figure(figsize=(16, 18))

    plt.subplot(2, 1, 1)
    for j in range(20):
        f_plot = []
        c_plot = []
        for f, c in zip(freq_indx, cl):
            if j < len(f):
                f_plot.append(f[j])
                c_plot.append(c[j])
        plt.plot(f_plot, c_plot, label=f'mode {j}')
    plt.legend()
    plt.xlabel('f [Hz]')
    plt.ylabel('c_L [m/s]')
    plt.title('figure 1: Love wave velocity as a function of frequency',
              fontsize=9)

    plt.subplot(2, 1, 2)
    for j in range(6):
        f_plot = []
        w_plot = []
        for f, w in zip(freq_indx, wvl):
            if j < len(f):
                f_plot.append(f[j])
                w_plot.append(w[j])
        plt.plot(f_plot, w_plot, label=f'mode {j}')
    plt.legend()
    plt.xlabel('f [Hz]')
    plt.ylabel('wavelength [m]')
    plt.title('figure 2: wavelength as a function of frequency', fontsize=9)

    plt.savefig("figures/mode")

if __name__ == "__main__":
    main()
