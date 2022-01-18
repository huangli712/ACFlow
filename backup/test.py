import scipy
import numpy
import scipy.linalg

A = numpy.array([[3.81631, 0.0, 1.6561],
                 [0.0, 7.9080, 0.0],
                 [1.6561, 0.0, 7.3933]], dtype = numpy.float64)
W, V = scipy.linalg.eigh(A)
#print(W)
print(V)