import numpy
import scipy.interpolate

mesh, image = numpy.loadtxt('Aout.data', unpack = True, usecols = (0,1))
simage = numpy.cumsum(image)
f = scipy.interpolate.interp1d(simage, mesh)

smax = numpy.amax(simage)
smin = numpy.amin(simage)
y = numpy.linspace(smin, smax, 1000)
x = f(y)

for i in range(1000):
    print(i, x[i])