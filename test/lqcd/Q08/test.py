import numpy
import scipy.interpolate
import sys

def time_to_freq_B(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
    """ forward fourier transformation from imaginary time to matsubara frequency,
        this function is only suitable for the bosonic function
        time_mesh : imaginary time mesh
        time_data : function on imaginary time axis
        freq_mesh : matsubara frequency mesh
        freq_data : function on matsubara frequency axis
        fast      : whether scipy.weave is used to accelerate the code
    """
    fit = scipy.interpolate.InterpolatedUnivariateSpline(time_mesh, time_data)
    ntime_dense = 4 * self.__ntime # used to build a dense imaginary time mesh
    # denser time mesh
    time_mesh_dense = Mesh.create_time_mesh(self.__beta, ntime_dense)
    # denser time data
    time_data_dense = fit(time_mesh_dense)
    for im in range(self.__nfreq):
        faux = time_data_dense * numpy.exp(1j * freq_mesh[im] * time_mesh_dense)
        # calculate \int^{\beta}_{0} f(\tau) e^{i \omega \tau} d\tau
        # now faux = f(\tau) e^{i \omega \tau}
        freq_data[im] = numpy.trapz(faux, time_mesh_dense)

print("hehe")
