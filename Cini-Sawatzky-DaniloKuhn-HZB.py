import numpy as np
from scipy.signal import fftconvolve, hilbert
from scipy.interpolate import interp1d
from scipy.optimize import brentq


class CSModel:
    """
    Cini-Sawatzky line shape model

    Ported from Igor Pro
    Danilo Kühn (HZB)
    Python conversion
    """

    def __init__(self, xdos, ydos):
        self.xdos = np.asarray(xdos)
        self.ydos = np.asarray(ydos)
        self.dx = self.xdos[1] - self.xdos[0]
        self.initialize()

    ###########################################################
    # Initialization
    ###########################################################

    def initialize(self):
        # normalize DOS
        self.ydos = self.ydos / np.sum(self.ydos * self.dx)

        # self convolution
        sc = fftconvolve(self.ydos,
                         self.ydos,
                         mode='same') * self.dx

        sc[sc < 1e-6] = 0

        # Hilbert transform
        analytic = hilbert(sc)
        self.ySCDOS = sc
        self.yHil = -np.pi * np.imag(analytic)

    ###########################################################
    # Gaussian
    ###########################################################
    @staticmethod
    def gaussian(x, mu, sigma):
        return np.exp(-(x-mu)**2/(2*sigma**2)) \
               /(np.sqrt(2*np.pi)*sigma)

    ###########################################################
    # Gaussian convolution
    ###########################################################
    def gaussconv(self, y, fwhm):
        if fwhm <= 0:
            return y.copy()

        sigma = fwhm / 2.355
        radius = int(np.ceil(6*sigma/self.dx))
        xx = np.arange(-radius, radius+1) * self.dx
        g = self.gaussian(xx, 0, sigma)
        g /= np.sum(g)

        return fftconvolve(y, g, mode="same")

    ###########################################################
    # Find bound state
    ###########################################################
    def bound_state(self, U):
        f = 1 - U*self.yHil
        s = np.sign(f)
        idx = np.where(np.diff(s))[0]

        if len(idx)==0:
            return None

        i = idx[-1]

        return brentq(
            lambda E: np.interp(E,
                                self.xdos,
                                f),
            self.xdos[i],
            self.xdos[i+1]
        )

    ###########################################################
    # Single CS peak
    ###########################################################
    def cs_peak(self,
                U,
                fwhm,
                intensity):
        helper = 1 - U*self.yHil
        y = self.ySCDOS / (
                helper**2
                + (np.pi*U*self.ySCDOS)**2
        )
        integral = np.sum(y)*self.dx
        y = self.gaussconv(y, fwhm)
        y *= intensity
        Eb = self.bound_state(U)

        if Eb is not None:
            sigma = fwhm/2.355
            weight = intensity*(1-integral)
            y += weight*self.gaussian(
                self.xdos,
                Eb,
                sigma
            )

        return y

    ###########################################################
    # Total model
    ###########################################################
    def spectrum(self,
                 x,
                 params):

        U1,U2,F1,F2,I1,I2,S1,S2 = params

        y1 = self.cs_peak(U1,F1,I1)
        y2 = self.cs_peak(U2,F2,I2)

        interp1 = interp1d(self.xdos + S1,
                           y1,
                           bounds_error=False,
                           fill_value=0)

        interp2 = interp1d(self.xdos + S2,
                           y2,
                           bounds_error=False,
                           fill_value=0)

        return interp1(x) + interp2(x)

# Example
import matplotlib.pyplot as plt

# Load DOS
xdos = np.loadtxt("dos_x.dat")
ydos = np.loadtxt("dos_y.dat")
model = CSModel(xdos, ydos)

x = np.linspace(0,25,2000)

params = [
    3.5,      # U1
    6.0,      # U2
    1.0,      # FWHM1
    1.3,      # FWHM2
    1.0,      # intensity1
    0.4,      # intensity2
    0.0,      # shift1
    3.5       # shift2
]

y = model.spectrum(x, params)

plt.plot(x,y)
plt.xlabel("Energy (eV)")
plt.ylabel("Intensity")
plt.show()
