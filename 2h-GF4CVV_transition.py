##################################################################
##################################################################
## ##
## Program for the calculation of the correlated two holes ##
## Green's function for a CVV Auger transition             ##
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
import datetime
from pprint import pprint
# return_pairs() creates all the possible combinations of the k_x
# and k_y values calcuted on the First Brillouin zone of the
# considered system or else of the CuO_2 planes.
def return_pairs(seq):
"""
return_pairs(seq) -> Given a sequence return a generator which
yields all the possible pairs, with repetitions.
Example:
>>> for item in return_pairs(xrange(2)): print item
(0, 0)
(0, 1)
(1, 0)
(1, 1)
"""
copy_of_seq = seq[:]
for i in seq:
for j in copy_of_seq:
yield i, j
# theta() is simply the definition by means of logic operators of
# the Heavyside function.
def theta(x):
""" theta (x) -> 1 if x < 0, 0 if x >= 0.
It is -1 times the Heavyside function.
"""
return x < 0
def my_module(x):
"""
my_modules (x) takes a value in R and brings it back into the
first Brillouin zone, -infinity < x < infinity -> - pi <= x < pi
"""
if scipy.sign(x) == 1:
while x > scipy.pi:
x = x-2*scipy.pi
else:
while x <= -1*scipy.pi:
x = x + 2*scipy.pi
return x
def my_range(start, stop, step):
"""
my_range(start, stop, step) yields values from start (included)
to stop incrementing by step. It is used for the omega vector:
that is, the energies where the Green's function will be evaluated.
It accepts any numeric values and not only integers like range and
xrange.
"""
value = start
if step <= 0:
print "Can only iterate forward! Step must be positive."
raise StopIteration
while value < stop:
yield value
value = value + step
# Main
class GreenFunction(object):
"""
__init__(precision, energies, hopping, coulomb_interactions,
fermi, label)
where precision is the number of digits requested,
energies is a dict with 'ed' and 'ep' keys,
hopping is a dict with 'tpd' and 'tpp' keys,
coulomb_interactions is a dict with 'udd' and 'upp' keys,
fermi is the fermi energy of the system.
tpd is the Cu-O hopping term;
tpp is the O-O hopping term;
ed is the 3d orbital binding energy
Since the oxygen 2p binding energy, ep, has been
set to zero, ed corresponds to the energy difference ed-ep
(we are in the hole picture);
label is a description of the system.
"""
def __init__(self, precision, energies, hopping,
coulomb_interactions, fermi, label):
self.precision = precision
# Parameter initialization
self.ed = energies['ed'] # d orbital
self.ep = energies['ep'] # p orbital
self.udd = coulomb_interactions['udd']
self.upp = coulomb_interactions['upp']
self.tpd = hopping['tpd']
self.tpp = hopping['tpp']
self.fermi = fermi
self.label = label
self.brillouin = [scipy.pi*(-1+2*n/L) for n in range(1,L+1)]
def hh(self, kx,ky):
# Definition of the system Hamiltonian, which has as variables
# the components of the wave vector k in the xOy plane and
# as parameters tpd, tpp and ed. Here:
ed = self.ed
ep = self.ep
tpd = self.tpd
tpp = self.tpp
return matrix([[ed,2*tpd*cos(kx/2),2*tpd*cos(ky/2)],
[2*tpd*cos(kx/2),ep, -4*tpp*cos(kx/2)*cos(ky/2)],
[2*tpd*cos(ky/2),-4*tpp*cos(kx/2)*cos(ky/2),ep]])
def compute(self, energy_range, interacting_green_function,
timeit,output=None):
"""
compute (energy_range, interacting_green_function, timeit,
output) -> It computes the correlated time ordered two
particles Green function, passed via interacting_green_function.
If timeit is True, it computes the time.
If output is an open file, it writes the results on it.
"""
self.kdata = self.eigensystem() # eigensystem for all the
# k vectors
final_data = []
print self.label
print "Job starting."
if timeit:
start = time.time() # Time counter
self.my_interacting = interacting_green_function
for omega in energy_range:
print '.',
qresult = self.compute_g0(omega)
final_data.append(self.compute_my_interacting(qresult))
if timeit:
end = time.time()
print "\nJob done in ", end - start, "seconds\n"
datapoints = [val.imag/(-1*math.pi*L**2) for val in final_data]
if output is None:
print datapoints
else:
output.write(self.label + '\n')
datapoints = ", ".join(str(value) for value in datapoints)
output.write(datapoints)
output.close()
def eigensystem(self):
"""
For each k vector in the brillouin zone, compute eigenvalues
and eigenvectors of the Hamiltonian.
Returns a dictionary.
"""
kdata={}
for index, kvector in enumerate(return_pairs(self.brillouin)):
kdata[index]=self.my_eig(*kvector)
return kdata
def my_eig(self, kx, ky):
"""
my_eig computes the eigenvalues and eigenvectors of the given
Hamiltonian and returns a Python dictionary:
kx, ky -> {0:{'eigval':value of the 1st eigenvalue,
'eigvec':1st eigenvector}, 1:{'eigval':value of the 2nd
eigenvalue, 'eigvec':2nd eigenvector}, 2:{'eigval':value of
the 3rd eigenvalue, 'eigvec':3rd eigenvector}} corresponding
to hh(kx, ky).
"""
precision = self.precision
res =eig(self.hh(kx,ky))
diz={}
for index in range(3):
diz[index]={'eigval':res[0][index],
'eigvec':scipy.around(scipy.asarray(res[1][:,index]),
decimals=precision)}
return diz
# compute_g0(omega) evaluates a 3*3 matrix for each considered
# (qx,qy) vector which represents the not correlated two holes
# Green's function of the system as a function of the hole
# binding energies omega.
def compute_g0(self, omega):
""" compute_g0(omega) evaluates a 3*3 matrix for each
(qx,qy) vector which represents the not correlated two
holes Green's function of the system as a function of the
hole binding energies omega.
"""
# To speed up lookup, bind to local names E_Fermi and scipy.outer
E_Fermi = self.fermi
out = scipy.outer
qsum = {}
for index, qvector in enumerate(return_pairs(self.brillouin)):
qx,qy = qvector
qsum[index]=0
for k_index, k_point in
enumerate(return_pairs(self.brillouin)):
kx, ky = k_point
qeig = self.my_eig(my_module(qx+kx),my_module(qy+ky))
for n in range(3):
for m in range(3):
v1=qeig[n]["eigvec"]*
self.kdata[k_index][m]["eigvec"]
num=out(v1,v1)*(1-theta(qeig[n]["eigval"]
-E_Fermi)-
theta(self.kdata[k_index][m]["eigval"]
-E_Fermi))
den=(omega - (qeig[n]["eigval"]
+self.kdata[k_index][m]["eigval"])
+1j*delta)*L**2
qsum[index] += num/den
return qsum
def compute_my_interacting(self, qresult):
"""
compute_GintCu() -> Calculates GintCu(a,b,c,d,m,f) for every
(qx,qy) in the first Brillouin zone and then sums all
these values
"""
Gnosum = 0
for qres in qresult.values():
Gnosum += self.my_interacting(qres[0,0],qres[0,1],
qres[0,2],qres[1,1],qres[1,2],qres[2,2], self.upp, self.udd)
return Gnosum
class ParameterSet(object):
pass
def cu_GreenInteracting(a,b,c,d,m,f, upp, udd):
"""
cu_GreenInteracting(a,b,c,d,m,f, upp, udd) ->
Returns the element [0,0] of the interacting
two holes Green's function matrix.
Its arguments are respectively:
a = g0[0,0], b = g0[0,1],c = g0[0,2],d = g0[1,1],m = g0[1,2],
f = g0[2,2]
"""
return (c*(c*upp - c*d*upp**2 + b*m*upp**2))/(1 - a*udd - d*upp
- f*upp - (b**2)*udd*upp - (c**2)*udd*upp + a*d*udd*upp
+ a*f*udd*upp + d*f*upp**2 - (m**2)*(upp**2)
+ (c**2)*d*udd*upp**2 + (b**2)*f*udd*upp**2
- a*d*f*udd*upp**2 - 2*b*c*m*udd*upp**2 + a*(m**2)*udd*upp**2)
+ (b*(b*upp - b*f*upp**2 + c*m*upp**2))/(1 -
a*udd - d*upp - f*upp - (b**2)*udd*upp - (c**2)*udd*upp
+ a*d*udd*upp + a*f*udd*upp + d*f*upp**2 - (m**2)*(upp**2)
+ (c**2)*d*udd*upp**2 + (b**2)*f*udd*(upp**2)
- a*d*f*udd*upp**2 - 2*b*c*m*udd*upp**2 +a*(m**2)*udd*upp**2)
+ (a*(1 - d*upp - f*upp + d*f*upp**2 - (m**2)*(upp**2)))
/(1 - a*udd - d*upp - f*upp - (b**2)*udd*upp -
(c**2)*udd*upp + a*d*udd*upp + a*f*udd*upp + d*f*upp**2
- (m**2)*(upp**2) + (c**2)*d*udd*upp**2 +(b**2)*f*udd*upp**2
- a*d*f*udd*upp**2 - 2*b*c*m*udd*upp**2 + a*(m**2)*udd*upp**2)
def o_GreenInteracting(a,b,c,d,m,f, upp, udd):
"""
o_GreenInteracting(a,b,c,d,m,f, upp, udd) ->
Returns the element [1,1] of the interacting
two holes Green's function matrix.
Its arguments are respectively:
a = g0[0,0], b = g0[0,1],c = g0[0,2],d = g0[1,1],m = g0[1,2],
f = g0[2,2]
"""
return (d*(1-a*udd-f*upp-(c**2)*udd*upp+a*f*udd*upp))/(1-a*udd-
d*upp-f*upp-(b**2)*udd*upp-(c**2)*udd*upp+a*d*udd*upp+
a*f*udd*upp+d*f*(upp**2)-(m**2)*(upp**2)+(c**2)*d*udd*(upp**2)+
(b**2)*f*udd*(upp**2)-a*d*f*udd*(upp**2)-2*b*c*m*udd*(upp**2)+
a*(m**2)*udd*(upp**2))+(m*(m*upp+b*c*udd*upp-a*m*udd*upp))/(1-
a*udd-d*upp-f*upp-(b**2)*udd*upp-(c**2)*udd*upp+a*d*udd*upp+
a*f*udd*upp+d*f*(upp**2)-(m**2)*(upp**2)+(c**2)*d*udd*(upp**2)+
(b**2)*f*udd*(upp**2)-a*d*f*udd*(upp**2)-2*b*c*m*udd*(upp**2)+
a*(m**2)*udd*(upp**2))+(b*(b*udd-b*f*udd*upp+c*m*udd*upp))/(1-
a*udd-d*upp-f*upp-(b**2)*udd*upp-(c**2)*udd*upp+a*d*udd*upp+
a*f*udd*upp+d*f*(upp**2)-(m**2)*(upp**2)+(c**2)*d*udd*(upp**2)+
(b**2)*f*udd*(upp**2)-a*d*f*udd*(upp**2)-2*b*c*m*udd*(upp**2)+
a*(m**2)*udd*(upp**2))
if __name__ == '__main__':
import psyco
psyco.full()
precision=5 # number of digits
L = 12 # L*L is the total number of k-vectors
delta = 0.1 # broadening of the calculated lineshape. It is
# equivalent to convolve the computed line shape
# with a Lorentzian curve having an HWHM of delta
green = {"Copper":cu_GreenInteracting, "Oxygen":o_GreenInteracting}
parameters = []
first_set = ParameterSet()
first_set.interacting = green["Copper"]
first_set.hopping = {'tpd':1.5, 'tpp': 0.6}
first_set.energies = {'ed': -3.3+(7.9/2)*0.273922,'ep':(3.6/2)*0.119984}
first_set.coulomb = {'udd':7.9,'upp': 3.6}
first_set.fermi = -5.2489
first_set.label = """
Description of the first parameter set.
"""
Description of the second parameter set.
"""
second_set.outputname = "today_Oxygen_HFBLA_L12_nh08472.dat"
parameters.append(first_set)
parameters.append(second_set)
for parameter in parameters:
nrg=my_range(-20.0,-10.0,0.5) # omega vector
my_greenfunction = GreenFunction(precision, parameter.energies,
parameter.hopping, parameter.coulomb, parameter.fermi,
parameter.label)
try:
today = str(datetime.datetime.today()).split()[0]
filename = today + '_' + parameter.outputname
my_file = open(filename, 'wb')
except AttributeError:
my_file = None
my_greenfunction.compute(nrg, parameter.interacting,
timeit=True, output=my_file)
