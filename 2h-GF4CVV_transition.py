##################################################################
##################################################################
## ##
## Program for the calculation of the correlated two holes ##
## Green's function for a CVV Auger transition  ##
## ##
##################################################################
##################################################################
# Everything is packed in classes to achieve better usability and
# performances.
# Import of the modules needed for the computations and the GUI
from __future__ import division
from numpy import matrix
from math import cos
import scipy
import numpy
try:
from pylab import axis, eig, plot, show, xlabel, ylabel, title
except:
print 'Couldn\'t find pylab, do something about that.'
from scipy.linalg import eig
import math
import time
