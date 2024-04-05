import numpy as np
import matplotlib.pyplot as plt
from goph420_lab02.root_finder import root_newton_raphson


def main():

    rho1 = 1800.0   # kg/m**3
    rho2 = 2500.0   # kg/m**3
    beta1 = 1900.0  # m/s
    beta2 = 3200.0  # m/s
    H = 4.0e03      # m
    undr_rt = H**2 * (beta1**-2-beta2**-2)

    zeta_max = np.sqrt(H**2 * (beta1**-2-beta2**-2))  # maximum value for zeta

    # initial root guess function
    def initial(frq, _k):
        return (((2*_k+1)/(4*frq))-1e-6)

    # list of frequencies in Hz to calculate initial guess at
    freq_list = np.array(np.linspace(1, 100, 50))
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
            return ((rho2/rho1)*np.sqrt(undr_rt - _zeta**2)/_zeta
                    - (np.tan(2*np.pi*freq*_zeta)))

        def dgdx(_zeta):
            return (-((2 * np.pi * freq) / np.cos(2 * freq * np.pi * _zeta)**2
                    + (rho2 / rho1) * np.sqrt(undr_rt - _zeta ** 2) / _zeta**2
                    + (rho2 / rho1) / (np.sqrt(undr_rt - _zeta ** 2))))

        if x0 > zeta_max:
            x0 = zeta_max - 1e-5

        while x0 < zeta_max:  # iterating over the modes to get to zeta max

            # Using the newton raphson root finding method
            zeta_k, _, _ = root_newton_raphson(x0, g, dgdx)

            # putting data of root for zeta's value
            # and the values that can then be derived from this value
            # into their respective lists
            zeta_mode.append(zeta_k)
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

    plt.figure(figsize=(10, 30))

    # plot frequency vs love wave velocity
    figplot(1, 21, freq_indx, cl, 'Frequency [Hz]', 'Love wave velocity [m/s]',
            'Figure 1: Love Wave Velocity as a Function of Frequency')
    # plot frequency vs wavelength
    figplot(2, 11, freq_indx, wvl,  'Frequency [Hz]', 'Wavelength [m]',
            'Figure 2: Wavelength as a Function of Frequency')
    # plot frequency vs zeta
    figplot(3, 11, freq_indx, wvl,   'Frequency [Hz]', 'Zeta',
            'Figure 3: Zeta as a Function of Frequency')

    plt.savefig("figures/mode")


def figplot(fignum, rng, x, y, x_lbl, y_lbl, title):

    plt.subplot(3, 1, fignum)
    for j in range(rng):
        x_plot = []
        y_plot = []
        for m, n in zip(x, y):
            if j < len(m):
                x_plot.append(m[j])
                y_plot.append(n[j])
        plt.plot(x_plot, y_plot, label=f'mode {j}')
    plt.legend()
    plt.xlabel(x_lbl)
    plt.ylabel(y_lbl)
    plt.title(title)


if __name__ == "__main__":
    main()
