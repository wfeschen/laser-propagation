import numpy as np
import matplotlib.pyplot as plt


def lorentzian(x, x0, gamma):
    """

    :param x: Parameter
    :param x0: Offset
    :param gamma: Width
    :return:
    """
    return (gamma**2 / ((x-x0)**2 + gamma**2))


def constant(x, const):
    return np.ones_like(x) * const


if __name__ == "__main__":
    delimiter = "   "
    z_start = 0.45
    z_end = 0.55
    Nz = 300
    lorentz_width = 100e-6
    p_max = 0.6

    nozzle_diameter = 500e-6
    nozzle_position = 5e-1
    nozzle_start = nozzle_position - nozzle_diameter/2.
    nozzle_end = nozzle_position + nozzle_diameter/2.

    z = np.linspace(z_start, z_end, Nz)
    pressure = np.zeros_like(z)
    z_before_nozzle = z[z < nozzle_start]
    z_nozzle = z[(z <= nozzle_end) * (nozzle_start <= z)]
    z_after_nozzle = z[z > nozzle_end]

    p_before_nozzle = lorentzian(z_before_nozzle, nozzle_start, lorentz_width)
    p_nozzle = constant(z_nozzle, 1)
    p_after_nozzle = lorentzian(z_after_nozzle, nozzle_end, lorentz_width)
    print(z_before_nozzle)
    print(z_nozzle)

    p = np.concatenate((p_before_nozzle, p_nozzle, p_after_nozzle)) * p_max

    np.savetxt("p.dat", np.column_stack([z, p]), delimiter=delimiter)

    plt.figure()
    # plt.plot(z_before_nozzle, p_before_nozzle)
    # plt.plot(z_nozzle, p_nozzle)
    # plt.plot(z_after_nozzle, p_after_nozzle)
    plt.plot(z, p)
    plt.show()
