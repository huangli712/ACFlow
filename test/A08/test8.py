import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, "/Users/lihuang/Working/working/ana_cont-master/")
import ana_cont.continuation as cont

true_w = np.linspace(-5.,5.,num=5001, endpoint=True)
beta = 40.
niw = 20
iw = np.pi/beta * (2.*np.arange(niw) + 1.)
def gauss_peak(center, width):
    return np.exp(-(true_w-center)**2/(2.*width**2)) / (width*np.sqrt(2.*np.pi))
true_spec_1 = 0.5*(gauss_peak(1., 0.2) + gauss_peak(2.,0.7))
true_spec_2 = 0.5*(gauss_peak(-1., 0.25) + gauss_peak(-2.1, 0.6))

spec_matrix = np.zeros((2,2,true_w.shape[0]))
spec_matrix[0,0] = true_spec_1
spec_matrix[1,1] = true_spec_2

rot_ang = 0.1
rot_mat = np.array([[np.cos(rot_ang), np.sin(rot_ang)],[-np.sin(rot_ang), np.cos(rot_ang)]])

true_spec_rot = np.einsum('ij,jkw,kl', rot_mat, spec_matrix, rot_mat.transpose())

#plt.plot(true_w, true_spec_1, color='black')
#plt.plot(true_w, true_spec_2, color='black')
#plt.plot(true_w, true_spec_rot[0,0], label='diag1', ls='--')
#plt.plot(true_w, true_spec_rot[1,1], label='diag2', ls='--')
#plt.plot(true_w, true_spec_rot[0,1], label='offd1', ls='--')
#plt.plot(true_w, true_spec_rot[1,0], label='offd2', ls='--')
#plt.legend()
#plt.show()

kernel = 1./(1j*iw[:,None] - true_w[None,:])
giw = np.trapz(kernel[None,None,:,:]*true_spec_rot[:,:,None,:], true_w, axis=3)
#plt.plot(iw, giw[1,1].real, ls='--')
#plt.plot(iw, giw[1,1].imag)
#plt.show()

wgrid = np.linspace(-4, 4, num=400)
#model = np.ones_like(wgrid)
#model = np.exp(-(wgrid-0.9)**2) * (wgrid-0.9)**2
model_diag = np.exp(-(wgrid/2.0)**2) / (2.0 * np.sqrt(np.pi))
model_diag /= np.trapz(model_diag, wgrid)
errfac = 0.00001
err = errfac * np.ones_like(iw)
#perfect_model = np.abs(np.interp(wgrid, true_w, true_spec))


f = open("green00.data", "w")
for i in range(niw):
    print(iw[i], " ", giw[0,0][i].real, giw[0,0][i].imag," ", 1e-5, file = f)
f.close()

f = open("green11.data", "w")
for i in range(niw):
    print(iw[i], " ", giw[1,1][i].real, giw[1,1][i].imag," ", 1e-5, file = f)
f.close()

f = open("green01.data", "w")
for i in range(niw):
    print(iw[i], " ", giw[0,1][i].real, giw[0,1][i].imag," ", 1e-5, file = f)
f.close()


probl_00 = cont.AnalyticContinuationProblem(im_axis=iw, re_axis=wgrid,
                                            im_data=giw[0,0], kernel_mode='freq_fermionic')
probl_11 = cont.AnalyticContinuationProblem(im_axis=iw, re_axis=wgrid,
                                            im_data=giw[1,1], kernel_mode='freq_fermionic')
probl_01 = cont.AnalyticContinuationProblem(im_axis=iw, re_axis=wgrid,
                                            im_data=giw[0,1], kernel_mode='freq_fermionic')

sol_00,sol_all_00 = probl_00.solve(method='maxent_svd', alpha_determination='chi2kink',
                                   optimizer='newton', stdev=err, model=model_diag,
                                   alpha_start=1e18, alpha_end=1e-10, fit_position=2.5,
                                   offdiag=False, interactive=True)
sol_11,sol_all_11 = probl_11.solve(method='maxent_svd', alpha_determination='chi2kink',
                                   optimizer='newton', stdev=err, model=model_diag,
                                   alpha_start=1e18, alpha_end=1e-10, fit_position=2.5,
                                   offdiag=False, interactive=True)

model_offd = np.sqrt(sol_00.A_opt * sol_11.A_opt)

f = open("model.data", "w")
for i in range(400):
    print(wgrid[i], " ", model_offd[i], file = f)
f.close()

#sys.exit(-1)


#plt.plot(true_w, true_spec_rot[0,0], color='blue', ls='--', label='true 1')
#plt.plot(true_w, true_spec_rot[1,1], color='red', ls='--', label='true 2')
#plt.plot(wgrid, sol_00.A_opt, color='blue', label='maxent 1')
#plt.plot(wgrid, sol_11.A_opt, color='red', label='maxent 2')
#plt.plot(wgrid, model_offd, color='green', label='default model for offdiag')
#plt.legend()
#plt.show()

#plt.plot(true_w, true_spec_rot[0,0], color='blue', ls='--')
#plt.plot(true_w, true_spec_rot[1,1], color='red', ls='--')
#plt.plot(wgrid, sol_00.A_opt, color='blue')
#plt.plot(wgrid, sol_11.A_opt, color='red')
#plt.plot(wgrid, model_offd, color='green')
#plt.xlim(-1,1)
#plt.show()

sol_01, sol_all_01 = probl_01.solve(method='maxent_svd', alpha_determination='chi2kink', optimizer='newton',
                                    stdev=err, model=model_offd, offdiag=True,
                                    preblur=False, blur_width=0.05, # preblur not necessary
                                    alpha_start=1e15, alpha_end=1e-15,
                                    interactive=True)

plt.plot(wgrid, sol_01.A_opt, color='green', label='maxent')
plt.plot(true_w, true_spec_rot[0,1], color='green', ls='--', label='true')
plt.legend()
plt.show()

plt.plot(wgrid, sol_01.A_opt, color='green', label='maxent')
plt.plot(true_w, true_spec_rot[0,1], color='green', ls='--', label='true')
plt.xlim(-1,1)
plt.legend()
plt.show()
