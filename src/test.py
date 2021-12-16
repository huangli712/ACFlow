import scipy
import numpy
import scipy.linalg

A = numpy.array([[3.81631842971968e-11, -0.0, 1.6561455252407786e-8],
                 [0.0, 7.907971584057321e-5, -0.0],
                 [1.6561455252407786e-8, 0.0, 7.393321992039924e-6]])
W, V = scipy.linalg.eigh(A)
print(W)
print(V)