import numpy as np
import matplotlib.pyplot as plt

PI = np.pi

# ---------------------------------------------------------
# numerical integration
# ---------------------------------------------------------
def rintegrate(f, dx):
    return np.sum(f) * dx

# ---------------------------------------------------------
# self convolution: 1H DOS --> 2H DOS
# ---------------------------------------------------------
def selfconvolute(d01, de):
    """
    d02(E) = integral d01(E') d01(E-E') dE'
    """
    d02 = np.convolve(d01, d01) * de
    return d02

# ---------------------------------------------------------
# Hilbert transform (simple version)
# ---------------------------------------------------------
def hilbert_transform(d, eta, ieglw, de, neg):
    """
    Equivalent to Fortran routine hilberttransf
    """
    ned = len(d)
    g = np.zeros(neg, dtype=complex)

    for ie in range(neg):
        je = np.arange(ned)
        denom = de * (ie + ieglw - je) - 1j * eta
        integrand = d / denom
        g[ie] = np.sum(integrand) * de

    return g


# ---------------------------------------------------------
# Cini-Sawatzky interacting Green function
# ---------------------------------------------------------
def cini_model(g02, U):
    return g02 / (1.0 - U * g02)


# =========================================================
# PARAMETERS
# =========================================================

ne = 100               # number of energy points
alpha = 1.0             # half bandwidth
bandwidth = 2 * alpha   

eta = 0.002   #! imaginary energy for evaluating Green functions
extraw = 2.0  #! extra plotting range outside the extrema of d01, d02
extra = int(extraw / (bandwidth / ne))

de = bandwidth / ne   #! energy increment from d01(1) to d01(2)
estart = -alpha       #! energy corresponding to d01(1)

# =========================================================
# 1-hole DOS (square band)
# =========================================================
d01 = np.ones(ne) / bandwidth

E1 = estart + np.arange(ne) * de

# =========================================================
# 2-hole DOS
# =========================================================
d02 = selfconvolute(d01, de)

E2 = 2 * estart + np.arange(len(d02)) * de

# =========================================================
# Green functions
# =========================================================
g01 = hilbert_transform(
    d01,
    eta,
    -extra,
    de,
    ne + 2 * extra
)

Eg1 = estart - extra * de + np.arange(len(g01)) * de


g02 = hilbert_transform(
    d02,
    eta,
    -extra,
    de,
    len(d02) + 2 * extra
)

Eg2 = 2 * estart - extra * de + np.arange(len(g02)) * de


# =========================================================
# Interacting spectrum
# =========================================================
U_values = [0.0, 1.0, 1.8, 3.0]

spectra = []

for U in U_values:
    g2 = cini_model(g02, U)
    spectrum = np.imag(g2) / PI
    spectra.append(spectrum)


# =========================================================
# PLOTS
# =========================================================

# ---------- input DOS ----------
plt.subplot(2, 2, 1) # 2 row, 2 columns, 1st subplot
plt.plot(E1, d01)
plt.xlabel("Energy")
plt.ylabel("1H DOS")
plt.title("Input one-hole DOS")
plt.grid()

# ---------- noninteracting 2-hole DOS ----------
plt.subplot(2, 2, 2) # 2 row, 2 columns, 2nd subplot
plt.plot(E2, d02)
plt.xlabel("Energy")
plt.ylabel("2H DOS")
plt.title("Non-interacting two-hole DOS")
plt.grid()

# ---------- Green function ----------
plt.subplot(2, 2, 3) # 2 row, 2 columns, 3rd subplot
plt.plot(Eg2, np.real(g02), label="Re G02")
plt.plot(Eg2, np.imag(g02), label="Im G02")
plt.xlabel("Energy")
plt.ylabel("G02")
plt.title("Non-interacting two-hole Green function")
plt.legend()
plt.grid()

# ---------- interacting spectra ----------
plt.subplot(2, 2, 4) # 2 row, 2 columns, 4th subplot

for U, spec in zip(U_values, spectra):
    plt.plot(Eg2, spec, label=f"U={U}")

plt.ylim(0, 1) 

plt.xlabel("Energy")
plt.ylabel("Im(G2)/π")
plt.title("Cini-Sawatzky spectra")
plt.legend()
plt.grid()

plt.show()
