import sys, os
import numpy as np

w_real = np.linspace(-5., 5., num=2001, endpoint=True)
spec_real = np.exp(-(w_real-0.5)**2 / (2.*0.2**2))
spec_real += 0.3 * np.exp(-(w_real+2.5)**2 / (2.*0.8**2))
spec_real /= np.trapz(spec_real, w_real) # normalization

beta = 10.
iw = np.pi/beta * (2.*np.arange(10) + 1.)

noise_amplitude = 1e-4 # create gaussian noise
rng = np.random.RandomState(123)
noise_abs = rng.normal(0., noise_amplitude, iw.shape[0])
noise_phase = rng.uniform(0., 2.*np.pi, iw.shape[0])
noise = noise_abs * np.exp(1j * noise_phase)
print(noise_abs)

kernel = 1./(1j*iw[:,None] - w_real[None,:])
gf_mats = np.trapz(kernel*spec_real[None,:], w_real, axis=1) + noise

err = np.ones_like(iw)*noise_amplitude

f = open("green.data1", "w")
for i in range(len(iw)):
    print(iw[i], " ", gf_mats[i].real, " ", gf_mats[i].imag, " ", err[i], file = f)
f.close()
