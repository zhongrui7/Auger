import numpy as np
from datetime import datetime

# Constants
PI = 3.14159265358979
r8 = np.float64  # Double precision equivalent

# Global variables
ne = 0  # Number of energy points in the band
d01 = None  # One-hole density of states (non-interacting)
d02 = None  # Two-hole density of states (non-interacting)
spec = None  # Spectrum: sum of all interacting spectra
g01 = None  # One-hole Green's function (non-interacting)
g02 = None  # Two-hole Green's function (non-interacting)
g2 = None   # Two-hole Green's function (interacting)
nelec = 0.0  # Number of electrons in the band
estart = 0.0  # Energy corresponding to d01[0]
de = 0.0      # Energy increment
extraw = 0.0  # Extra plotting range
eta = 0.0     # Imaginary energy for Green's functions
interaction = 0.0  # Hole-hole interaction strength
augerm = 0.0  # Auger parameter (?)
norm = 0.0    # Norm of density of states
bandwidth = 0.0  # Width of the one-hole band
alpha = 0.0      # Parameter alpha (bandwidth/2)
gamma = 0.0      # Parameter gamma (W/alpha)
extra = 0        # Number of extra energy points
type_ = ""       # Calculation type: 'file' or 'square'
inpfile = ""     # Input file for 1P-DOS

def selfconvolute(d01, de):
    """
    Computes the self-convolution of d01 to produce d02.
    Requires len(d02) = 2 * len(d01) + 1.
    """
    ne = len(d01)
    d02 = np.zeros(2 * ne + 1, dtype=r8)
    for ie in range(2 * ne + 1):
        jemin = max(0, ie - ne)
        jemax = min(ne, ie)
        aux = np.zeros(jemax - jemin, dtype=r8)
        for je in range(jemin, jemax):
            aux[je - jemin] = d01[je] * d01[ie - je]
        d02[ie] = rintegrate(aux, de)
    return d02

def hilberttransf(d, eta, ieglw, de):
    """
    Computes the Hilbert transform of the density of states d to produce Green's function g.
    """
    ned = len(d)
    neg = len(g01) if 'g01' in globals() else len(g02)
    g = np.zeros(neg, dtype=np.complex128)
    for ie in range(neg):
        integrand = np.zeros(ned, dtype=np.complex128)
        for je in range(ned):
            integrand[je] = d[je] / (de * (ie + ieglw - je) - 1j * eta)
        g[ie] = cintegrate(integrand, de)
    return g

def cinimodel(g02, interaction):
    """
    Computes the interacting two-hole Green's function from the non-interacting one.
    """
    return g02 / (1.0 - interaction * g02)

def rintegrate(funct, dx):
    """
    Integrates a real-valued function using simple summation.
    """
    return np.sum(funct) * dx

def cintegrate(funct, dx):
    """
    Integrates a complex-valued function using simple summation.
    """
    return np.sum(funct) * dx

def rplot(vect, start, de, un):
    """
    Writes real-valued vector to file with energies.
    """
    with open(f"output_{un}.dat", 'w') as f:
        for i, val in enumerate(vect):
            f.write(f"{(i) * de + start:15.7f} {val:15.7e}\n")

def cplot(vect, start, de, un):
    """
    Writes complex-valued vector to file with energies.
    """
    with open(f"output_{un}.dat", 'w') as f:
        for i, val in enumerate(vect):
            f.write(f"{(i) * de + start:15.7f} {val.real:15.7e} {val.imag:15.7e}\n")

def main():
    global ne, d01, d02, spec, g01, g02, g2, nelec, estart, de, extraw, eta
    global interaction, augerm, norm, bandwidth, alpha, gamma, extra, type_, inpfile

    # Input calculation type
    while type_ not in ['file', 'square']:
        type_ = input("DOS type [ file | square ]: ").strip()

    # Input number of energy points
    ne = int(input("number of energy points in the band: "))

    # Allocate d01
    d01 = np.zeros(ne, dtype=r8)

    if type_ == 'square':
        alpha = float(input("alpha (bandwidth/2): "))
        bandwidth = 2 * alpha
        d01[:] = 1.0 / bandwidth
        de = bandwidth / ne
        estart = -alpha
    else:  # type_ == 'file'
        inpfile = input("file name for 1P-DOS: ").strip()
        d01 = np.loadtxt(inpfile)[::-1]  # Reverse DOS
        nelec = float(input("number of electrons: "))
        d01 = d01 / nelec
        estart = float(input("energy start: "))
        de = float(input("energy increment: "))
        estart = -(estart + de * ne)

    extraw = float(input("extra plotting range: "))
    extra = int(extraw / de)
    eta = float(input("imaginary energy: "))

    # Allocate arrays
    g01 = np.zeros(ne + 2 * extra, dtype=np.complex128)
    d02 = np.zeros(2 * ne + 1, dtype=r8)
    g02 = np.zeros(2 * ne + 1 + 2 * extra, dtype=np.complex128)
    g2 = np.zeros(2 * ne + 1 + 2 * extra, dtype=np.complex128)
    spec = np.zeros(2 * ne + 1 + 2 * extra, dtype=r8)

    # Welcome message
    print("\n" + "="*80 + "\n")
    print("     Welcome to Cini-Sawatzky Model for CVV Auger lineshape...")
    print("\n")

    # Non-interacting 1H-DOS
    rplot(d01, estart, de, 10)
    norm = rintegrate(d01, de)
    print(f"1H-DOS norm: {norm:20.12e}")

    # Non-interacting 2H-DOS
    d02 = selfconvolute(d01, de)
    rplot(d02, 2 * estart, de, 20)
    norm = rintegrate(d02, de)
    print(f"2H-DOS norm: {norm:20.12e}")

    # Non-interacting 1H-Green function
    g01 = hilberttransf(d01, eta, -extra, de)
    cplot(g01, estart - extra * de, de, 11)
    norm = rintegrate(g01.imag, de)
    print(f"1H-DOS norm (from G^0_1): {norm / PI:20.12e}")

    # Non-interacting 2H-Green function
    g02 = hilberttransf(d02, eta, -extra, de)
    cplot(g02, 2 * estart - extra * de, de, 21)
    norm = rintegrate(g02.imag, de)
    print(f"2H-DOS norm (from G^0_2): {norm / PI:20.12e}")

    # Interacting 2H-Green function
    spec[:] = 0.0
    for i in range(1000):
        augerm = 1.0
        if type_ == 'square':
            gamma = float(input("h-h: gamma (=W/alpha) [<0 to exit]: "))
            interaction = alpha * gamma
        else:
            interaction = float(input("h-h: interaction [<0 to exit]: "))
            augerm = float(input("augerm (?): "))
        if interaction < 0:
            break
        g2 = cinimodel(g02, interaction)
        rplot(g2.imag / PI, 2 * estart - extra * de, de, 3000 + i)
        norm = rintegrate(g2.imag, de)
        print(f"{i:3d}: interacting 2H-DOS norm: {norm / PI:20.12e}")
        spec += augerm * g2.imag / PI

    rplot(spec, 2 * estart - extra * de, de, 40)

    print("\n" + "="*76 + " Bye\n")

if __name__ == "__main__":
    main()
