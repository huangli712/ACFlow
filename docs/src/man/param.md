*A comprehensive dictionary about parameters.*

## [Overview](@id param)

The official configuration file for the ACFlow toolkit is `case.toml`. This page contains all the valid parameters that can appear in `case.toml`. As for the format of `case.toml`, please look at [`case.toml`](input.md).

```@contents
Pages = ["param.md"]
Depth = 2:3
```

## [[BASE] Block](@id base_block)

!!! note

    This block is mandatory. The parameters in this block are useful for all the solvers.

### [finput](@id finput)

*Definition:*

> Filename for the input data. The input data should be stored in a column-wised and formated (CSV-like) text file.

*Type:*

> String.

*Example:*

> finput = "gtau.data"

*Comment:*

> This parameter is mandatory.

### [solver](@id solver)

*Definition:*

> This parameter specifies the solvers that used to solve the analytic continuation problem. Now the ACFlow toolkit supports seven different solvers. They are as follows:
>
> * MaxEnt
> * BarRat
> * NevanAC
> * StochAC
> * StochSK
> * StochOM
> * StochPX
>
> Here, `MaxEnt` means the maximum entropy method. The `MaxEnt` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "MaxEnt"`, then the `[MaxEnt]` block must be available in the configuration file.
>
> `BarRat` means the barycentric rational function approximation. The `BarRat` solver can be used to treat the correlators in Matsubara frequency axis only. If `solver = "BarRat"`, then the `[BarRat]` block must be available in the configuration file.
>
> `NevanAC` means the Nevanlinna analytical continuation. The `NevanAC` solver can be used to treat the fermionic correlators in Matsubara frequency. Note that this solver is extremely sensitive to the noise level of the input data.
>
> `StochAC` means the stochastic analytic continuation method (K. S. D. Beach's algorithm). The `StochAC` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochAC"`, then the `[StochAC]` block must be available in the configuration file.
>
> `StochSK` means the stochastic analytic continuation method (A. W. Sandvik's algorithm). The `StochSK` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochSK"`, then the `[StochSK]` block must be available in the configuration file.
>
> `StochOM` means the stochastic optimization method. The `StochOM` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochOM"`, then the `[StochOM]` block must be available in the configuration file.
>
> `StochPX` means the stochastic pole expansion method. The `StochPX` solver can be used to treat the correlators in Matsubara frequency axis only. If `solver = "StochPX"`, then the `[StochPX]` block must be available in the configuration file.

*Type:*

> String.

*Example:*

> solver = "MaxEnt"

*Comment:*

> This parameter is mandatory. Overall, we recommend the `MaxEnt`, `BarRat`, and `StochPX` solvers.

### [ktype](@id ktype)

*Definition:*

> It denotes the type of kernel functions. Now the ACFlow toolkit supports three types of kernel functions. They are:
>
> * fermi
> * boson
> * bsymm
>
> Here, `fermi` means fermionic kernel function, which reads
>
> ```math
> K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}},
> ```
>
> and
>
> ```math
> K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega}.
> ```
>
> `boson` means bosonic kernel function, which reads
>
> ```math
> K(\tau,\omega) = \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}},
> ```
>
> and
>
> ```math
> K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega}.
> ```
>
> `bsymm` means symmetric bosonic kernel function, which reads
>
> ```math
> K(\tau,\omega) = \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]} {1 - e^{-\beta\omega}},
> ```
>
> and
>
> ```math
> K(\omega_n, \omega) = \frac{-2\omega^2}{\omega_n^2 + \omega^2}.
> ```
>
> As for detailed formula for these kernel functions, please refer to the comments in `src/kernel.jl`.

*Type:*

> String.

*Example:*

> ktype = "fermi"

*Comment:*

> This parameter is mandatory. It must be compatible with the `grid` parameter.

### [mtype](@id mtype)

*Definition:*

> It denotes the type of default model functions. Now the ACFlow toolkit supports the following choices:
>
> * flat
> * gauss
> * 1gauss
> * 2gauss
> * lorentz
> * 1lorentz
> * 2lorentz
> * risedecay
> * file
>
> Here, `flat` means the flat model (i.e., constant), `gauss` means the Gaussian model, `1gauss` means the Shifted Gaussian model, `2gauss` means the Two Gaussians model, `lorentz` means the Lorentzian model, `1lorentz` means the Shifted Lorentzian model, `2lorentz` means the Two Lorentzians model, and `risedecay` means the Rise-And-Decay model.
>
> Besides `flat` and `file`, all the other model functions need additional parameters to customize them (Of course, the ACFlow toolkit will supplement default parameters). The parameters can be specified by the [`pmodel`](@ref pmodel) parameter.
>
> Especially, if `mtype = "file"`, then the default model function is encoded in `model.inp`. ACFlow will read this file and initialize the default model function automatically. Be careful, the mesh for this model function must be consistent with the one used in the analytic continuation calculations.
>
> As for detailed formula for these models, please refer to the comments in `src/model.jl`.

*Type:*

> String.

*Example:*

> mtype = "flat"

*Comment:*

> This parameter is mandatory. Only the `MaxEnt` solver need these model functions. The `StochAC` solver only supports the `flat` model. The `StochSK`, `StochOM`, and `StochPX` solvers are free of model functions.

### [grid](@id grid)

*Definition:*

> This parameter specifies the grid's type for input data in imaginary axis. Now the ACFlow toolkit supports the following choices:
>
> * ftime
> * fpart
> * btime
> * bpart
> * ffreq
> * ffrag
> * bfreq
> * bfrag
>
> Here, `ftime` means fermionic and imaginary time, `btime` means bosonic and imaginary time, `ffreq` means fermionic and Matsubara frequency, and `bfreq` means bosonic and Matsubara frequency. `fpart` means fermionic and imaginary time as well, but the grid of imaginary time might be incomplete. `bpart` is similar to `fpart`, but it is for the bosonic case. `ffrag` means fermionic and Matsubara frequency as well, but the grid of Matsubara frequency might be incomplete. `bfrag` is similar to `ffrag`, but it is for the bosonic case.

*Type:*

> String.

*Example:*

> grid = "ftime"

*Comment:*

> This parameter is mandatory. It must be compatible with the `ktype` parameter. If grid is "bfrag", the first Matsubara frequency point, i.e. ``i\omega_0 = 0``, should be kept. See also [`ngrid`](@ref ngrid).

!!! warning

    If the `BarRat` solver is employed, the `grid` parameter should be "ffreq", "ffrag", "bfreq", or "bfrag".

!!! warning

    If the `NevanAC` solver is employed, the `grid` parameter should be "ffreq" or "ffrag".

!!! warning

    If the `StochPX` solver is employed, the `grid` parameter should be "ffreq", "ffrag", "bfreq", or "bfrag".

### [mesh](@id mesh)

*Definition:*

> This parameter specifies the mesh's type for output data (usually the spectral functions) in real axis. Now the ACFlow toolkit supports the following choices:
>
> * linear
> * tangent
> * lorentz
> * halflorentz
>
> Here, `linear` means the Linear mesh, `tangent` means the Tangent mesh, `lorentz` means the Lorentzian mesh, and `halflorentz` means the Half-Lorentzian mesh.
>
> Notes that only the `linear` mesh is uniform, the other three meshes are non-uniform. And the `halflorentz` mesh is defined on the positive-half axis only.

*Type:*

> String.

*Example:*

> mesh = "linear"

*Comment:*

> This parameter is mandatory. See also [`nmesh`](@ref nmesh).

!!! warning

    As for the `StochOM` solver , it seems that the linear mesh works better.

### [ngrid](@id ngrid)

*Definition:*

> Number of grid points. The parameter, together with the `beta` and `grid` parameters, controls the generation of grid for input data.

*Type:*

> Integer.

*Example:*

> ngrid = 10

*Comment:*

> This parameter is mandatory. It must be compatible with the input data. See also [`grid`](@ref grid).

!!! warning

    If the `BarRat` solver is employed, the `ngrid` parameter should be large (~100-200).

### [nmesh](@id nmesh)

*Definition:*

> Number of mesh points. The parameter, together with the `wmax`, `wmin`, and `mesh` parameters, controls the generation of mesh for output data.

*Type:*

> Integer.

*Example:*

> nmesh = 501

*Comment:*

> This parameter is mandatory. See also [`mesh`](@ref mesh).

### [wmax](@id wmax)

*Definition:*

> Right boundary (maximum value) of mesh. Note that `wmax` should be always greater than `wmin`.

*Type:*

> Float.

*Example:*

> wmax = 10.0

*Comment:*

> This parameter is mandatory.

### [wmin](@id wmin)

*Definition:*

> Left boundary (minimum value) of mesh. Note that `wmax` should be always greater than `wmin`.

*Type:*

> Float.

*Example:*

> wmin = -10.0

*Comment:*

> This parameter is mandatory.

!!! warning

    If the `ktype = "bsymm"`, the `wmin` parameter should be 0.0. In other words, the spectral density is defined on the half positive axis.

### [beta](@id beta)

*Definition:*

> Inverse temperature ``\beta``. It is equal to ``1/T``.

*Type:*

> Float.

*Example:*

> beta = 10.0

*Comment:*

> This parameter is mandatory. This parameter must be compatible with the input data and grid. Specifically, for the imaginary time axis, the last grid point should be ``\beta``. As for the Matsubara frequency axis, the difference between two successive grid points should be ``\pi/\beta``.

### [offdiag](@id offdiag)

*Definition:*

> Is the input correlator the offdiagonal part in matrix-valued function? As for the offdiagonal correlator, the corresponding spectral function might be not positive-definite. Some tricks have been implemented to cure this issue.

*Type:*

> Bool.

*Example:*

> offdiag = false

*Comment:*

> This parameter is mandatory. This parameter is useful for the `MaxEnt` and `StochPX` solvers only.

!!! warning

    Now only the `MaxEnt` and `StochPX` solvers supports this parameter. On the other hand, the `auxiliary Green's function` algorithm works always for the solvers that don't support this parameter. The `BarRat` solver support analytic continuations for off-diagonal Green's functions, but it will ignore this parameter.

### [fwrite](@id fwrite)

*Definition:*

> Are the analytic continuation results written into external files? If it is false, then only the terminal output is retained and all the other outputs are disable. By default (if this parameter is missing or true), the files should be generated.

*Type:*

> Bool.

*Example:*

> fwrite = false

*Comment:*

> This parameter is optional.

### [pmodel](@id pmodel)

*Definition:*

> Additional parameters for customizing the model functions. Note that the `gauss`, `lorentz`, and `risedecay` models need one parameter ``\Gamma``. The `1gauss` and `1lorentz` models need two parameters, ``\Gamma`` and ``s``. The `2gauss` and `2lorentz` models need three parameters, ``\Gamma``, ``s_1``, and ``s_2``.
>
> The `pmodel` parameter is used to define these parameters. If there is only one element in `pmodel`, then ``\Gamma`` = `pmodel[1]`. If there are two elements in `pmodel`, then ``\Gamma`` = `pmodel[1]` and ``s`` = `pmodel[2]`. If there are three elements in `pmodel`, then ``\Gamma`` = `pmodel[1]`, ``s_1`` = `pmodel[2]`, and ``s_2`` = `pmodel[3]`.

*Type:*

> Array.

*Example:*

> pmodel = [1.0]

*Comment:*

> This parameter is optional.  The default values for ``\Gamma``, ``s``, ``s_1``, and ``s_2`` are 2.0, 2.0, -2.0, and 2.0, respectively.

### [pmesh](@id pmesh)

*Definition:*

> Additional parameters for customizing the mesh. The `tangent` mesh needs the ``f_1`` parameter. The `lorentz` and `halflorentz` meshes need the `cut` parameter. The `pmesh` parameter can be used to setup the two parameters. If `pmesh` contains one element or more than one elements, then ``f_1 \equiv `` `cut` ``\equiv`` `pmesh[1]`.

*Type:*

> Array.

*Example:*

> pmesh = [2.1]

*Comment:*

> This parameter is optional. The default values for ``f_1`` and `cut` are 2.1 and 0.01, respectively. See also [`mesh`](@ref mesh).

### [exclude](@id exclude)

*Definition:*

> Restriction of the energy range of the calculated spectral functions. This features is implemented by the `StochAC`, `StochSK`, `StochOM`, and `StochPX` solvers. In these solvers, the ``\delta`` or `box` functions, which are used to mimic the spectral functions, are restricted to live out of the given energy ranges. For example, `exclude = [[8.0,16.0]]` means that the energy range `[8.0,16.0]` is strictly forbidden.

*Type:*

> Array.

*Example:*

> exclude = [[-8.0,-4.0],[4.0,8.0]]

*Comment:*

> This parameter is optional. If you are using the `MaxEnt` solver, this parameter will be ignored. If solver = `StochPX` and offdiag = true, this parameter is mandatory. In this case, it is used to restrict the regions that the poles with positive weights can survive (or equivalently, the regions that the poles with negative weights can survice are also determined). For example, if exclude = [[-3.0,3.0]], wmin = -5.0, and wmax = 5.0, then the regions for poles with negative weights are [-3.0,3.0], while the regions for poles with positive weights are [-5.0,-3.0] U [3.0,5.0].

## [[MaxEnt] Block](@id maxent_block)

!!! note

    The parameters in this block are valid for the `MaxEnt` solver only.

!!! warning

    If `solver = "MaxEnt"`, the `[MaxEnt]` block must be available.

### [method](@id maxent_method)

*Definition:*

> How to determine the optimized ``\alpha`` parameter? The `MaxEnt` solver supports four different algorithms. They are
>
> * historic
> * classic
> * bryan
> * chi2kink
>
> Usually, the `chi2kink` algorithm is preferred.

*Type:*

> String.

*Example:*

> method = "bryan"

*Comment:*

> This parameter is mandatory. As for the underlying principles of these algorithms, please see [Maximum Entropy Method](@ref mem).

### [stype](@id maxent_stype)

*Definition:*

> Type of the entropic factor. The `MaxEnt` solver supports two schemes. They are
>
> * sj
> * br
>
> Here, `sj` means the Shannon-Jaynes entropy, while `br` means the Bayesian Reconstruction entropy. Usually, the Shannon-Jaynes entropy is preferred, since with it the positivity of the generated spectrum is always guaranteed. The Bayesian Reconstruction entropy tends to yield sharp features.

*Type:*

> String.

*Example:*

> stype = "sj"

*Comment:*

> This parameter is mandatory. As for the underlying principles of these entropic factors, please see [Maximum Entropy Method](@ref mem).

### [nalph](@id maxent_nalph)

*Definition:*

> Total number of the chosen ``\alpha`` parameters.

*Type:*

> Integer.

*Example:*

> nalph = 12

*Comment:*

> This parameter is mandatory. Only the `chi2kink` algorithm needs this parameter to control the number of ``\alpha`` parameters.

### [alpha](@id maxent_alpha)

*Definition:*

> Starting value for the ``\alpha`` parameter. The `MaxEnt` solver always starts with a huge ``\alpha`` parameter, and then decreases it gradually.

*Type:*

> Float.

*Example:*

> alpha = 1e9

*Comment:*

> This parameter is mandatory. It should be a very large number, such as ``10^9 \sim 10^{13}``.

### [ratio](@id maxent_ratio)

*Definition:*

> Scaling factor for the ``\alpha`` parameter. The next ``\alpha`` is equal to the current ``\alpha`` divided by `ratio`.

*Type:*

> Float.

*Example:*

>ratio = 10.0

*Comment:*

> This parameter is mandatory. It muse be larger than 1.0.

### [blur](@id maxent_blur)

*Definition:*

> Sometimes, the kernel functions and spectral functions can be preblurred to obtain smoother results. Shall we preblur them? If `blur` is larger than zero, then it means the blur parameter. If `blur` is smaller than zero, then it means that the preblur feature is disable.

*Type:*

> Float.

*Example:*

> blur = -1.0

*Comment:*

> This parameter is mandatory.

## [[BarRat] Block](@id barrat_block)

!!! note

    The parameters in this block are valid for the `BarRat` solver only.

!!! warning

    If `solver = "BarRat"`, the `[BarRat]` block must be available.

!!! warning

    The `BarRat` solver is still in development. Please use it at your own risk.

### [atype](@id barrat_atype)

*Definition:*

> Possible type of the spectrum.
>
> * cont
> * delta
>
> If it is `cont`, it means that the spectrum should be board and continuous. If it is `delta`, it means that the spectrum consists a few ``\delta``-like peaks. The `BarRat` solver will deduce the positions of the poles from the barycentric rational function approximation, and then the BFGS algorithm is used to determine the weights / amplitudes of these poles. The original and real-frequency Green's function are then reconstructed by using the pole representation.

*Type:*

> String.

*Example:*

> atype = "cont"

*Comment:*

> This parameter is mandatory. If `atype` is "delta", then the `pcut` and `eta` parameters will take effect. On the contrary, if `atype` is "cont", then the `pcut` and `eta` parameters will be ignored.

### [denoise](@id barrat_denoise)

*Definition:*

> This parameter specifies how to denoise the input data.
>
> * none
> * prony\_s
> * prony\_o
>
> The BarRat solver will adopt the Prony approximation to approximate the Matsubara data and suppress the noise. The `denoise` parameter is used to control whether the Prony approximation is actived. If it is "none", the Prony approximation is disabled. If it is "prony\_s", the Prony approximation will run once, and its accuracy is controlled by the `epsilon` parameter. If it is "prony\_o", an optimal Prony approximation is automatically determined. In such a case, the `epsilon` parameter is nonsense.
>
> If the Prony approximation is activated, the grid parameter should not be `ffrag` or `bfrag`.

*Type:*

> String.

*Example:*

> denoise = "none"

*Comment:*

> This parameter is mandatory. If the noise level is obvious, please set `denoise` to "prony\_s" or "prony\_o". If the noise level is small (such as ``\epsilon < 10^{-8}``), Prony approximation could lead to worse results.

### [epsilon](@id barrat_epsilon)

*Definition:*

> Threshold for the Prony approximation. It is used to control the accuracy of the Prony approximation. It can be considered as a measurement of the noise level of the input Matsubara data.

*Type:*

> Float.

*Example:*

> epsilon = 1e-10

*Comment:*

> This parameter is mandatory. But it is only useful when `denoise` is not "none". See [`denoise`](@ref barrat_denoise) for more details.
>
> `epsilon` should be set to the noise level of the input Matsubara data. It should not be too small or too large. In principle, `` \sigma[1] < \textrm{epsilon} < \sigma[\textrm{end}]``, where ``\sigma`` are the singular values.

### [pcut](@id barrat_pcut)

*Definition:*

> Cutoff for unphysical poles. This parameter is used to filter the unphysical poles generated by the AAA algorithm. Given a pole, if the imaginary part of its location or weight is larger than `pcut`, it should be discarded or removed.
>
> Sometimes if `pcut` is too small, all of the poles are consiered as unphysical and removed. The ACFlow will throw an error and stop. To cure this problem, please increase `pcut` and redo the calculation.

*Type:*

> Float.

*Example:*

> pcut = 1e-3

*Comment:*

> This parameter is mandatory. But it is only useful when `atype = "delta"`. See [`atype`](@ref barrat_atype) for more details.

### [eta](@id barrat_eta)

*Definition:*

> Tiny distance from the real axis. It is used to construct the retarded Matsubara Green's function within the pole representation.
>
> If eta is smaller than 1.0, then
>
> ```math
> G(\omega) = \sum_j \frac{\operatorname{Re} A_j}{\omega - \operatorname{Re} P_j + i\eta }
> ```
>
> If eta is larger than 1.0, then
>
> ```math
> G(\omega) = \sum_j \frac{A_j}{\omega - \operatorname{Re} P_j + i(\eta - 1) }
> ```

*Type:*

> Float.

*Example:*

> eta = 1e-2

*Comment:*

> This parameter is mandatory. But it is only useful when `atype = "delta"`. See [`atype`](@ref barrat_atype) for more details.
>
> Well, for normal Green's function, the imaginary parts of `A` and `P` must be quite small. So please let `eta < 1.0`. But for anormal Green's function, the imaginary parts of `A` and `P` might be quite large. You can set `eta > 1.0`. Anyway, more tests are essential.

## [[NevanAC] Block](@id nevanac_block)

!!! note

    The parameters in this block are valid for the `NevanAC` solver only.

!!! warning

    If `solver = "NevanAC"`, the `[NevanAC]` block must be available.

!!! warning

    This solver is numerically unstable, so use it at your own risk.

### [pick](@id nevanac_pick)

*Definition:*

> Check the Pick criterion or not. If `pick` is true, ACFlow will try to figure out the optimal number of the input data (i.e., how many data points are retained for further postprocessing) by using the Pick criterion.

*Type:*

> Bool.

*Example:*

> pick = true

*Comment:*

> This parameter is mandatory.

### [hardy](@id nevanac_hardy)

*Definition:*

> Perform Hardy basis optimization or not. The spectrum obtained by the Nevanlinna analytical continuation is usually wiggly. So, the Hardy basis optimization can help us smooth the spectrum.

*Type:*

> Bool.

*Example:*

> hardy = true

*Comment:*

> This parameter is mandatory. See also [`hmax`](@ref nevanac_hmax).

### [hmax](@id nevanac_hmax)

*Definition:*

> Upper cut off of Hardy order. In principle, the larger the Hardy order is, the smoother the obtained spectrum is. Usually `hmax = 20` is enough to get smooth spectrum.

*Type:*

> Integer.

*Example:*

> hmax = 50

*Comment:*

> This parameter is mandatory. See also [`hardy`](@ref nevanac_hardy).

### [alpha](@id nevanac_alpha)

*Definition:*

> Regulation parameter for smooth norm. If the `alpha` parameter is too large, the detailed features in the spectra could be smeared out.

*Type:*

> Float.

*Example:*

> alpha = 1e-4

*Comment:*

> This parameter is mandatory.

### [eta](@id nevanac_eta)

*Definition:*

> Tiny distance from the real axis. It is used to construct the Hardy matrix, instead of the Green's function.

*Type:*

> Float.

*Example:*

> eta = 1e-2

*Comment:*

> This parameter is mandatory.

## [[StochAC] Block](@id stochac_block)

!!! note

    The parameters in this block are valid for the `StochAC` solver only.

!!! warning

    If `solver = "StochAC"`, the `[StochAC]` block must be available.

### [nfine](@id stochac_nfine)

*Definition:*

> Number of points of a very fine linear mesh. This mesh is for the ``\delta`` functions.

*Type:*

> Integer.

*Example:*

> nfine = 10000

*Comment:*

> This parameter is mandatory.

### [ngamm](@id stochac_ngamm)

*Definition:*

> Number of ``\delta`` functions. Their superposition is used to mimic the spectral functions.

*Type:*

> Integer.

*Example:*

> ngamm = 512

*Comment:*

> This parameter is mandatory.

### [nwarm](@id stochac_nwarm)

*Definition:*

> Number of Monte Carlo thermalization steps.

*Type:*

> Integer.

*Example:*

> nwarm = 4000

*Comment:*

> This parameter is mandatory.

### [nstep](@id stochac_nstep)

*Definition:*

> Number of Monte Carlo sweeping steps.

*Type:*

> Integer.

*Example:*

> nstep = 4000000

*Comment:*

> This parameter is mandatory.

### [ndump](@id stochac_ndump)

*Definition:*

> Intervals for monitoring Monte Carlo sweeps. For every `ndump` steps, the `StochAC` solver will try to output some useful information to help diagnosis.

*Type:*

> Integer.

*Example:*

> ndump = 40000

*Comment:*

> This parameter is mandatory.

### [nalph](@id stochac_nalph)

*Definition:*

> Total number of the chosen ``\alpha`` parameters.

*Type:*

> Integer.

*Example:*

> nalph = 20

*Comment:*

> This parameter is mandatory.

### [alpha](@id stochac_alpha)

*Definition:*

> Starting value for the ``\alpha`` parameter. The `StochAC` solver always starts with a small ``\alpha`` parameter, and then increases it gradually.

*Type:*

> Float.

*Example:*

> alpha = 1.0

*Comment:*

> This parameter is mandatory.

### [ratio](@id stochac_ratio)

*Definition:*

> Scaling factor for the ``\alpha`` parameter. It should be larger than 1.0.

*Type:*

> Float.

*Example:*

> ratio = 1.2

*Comment:*

> This parameter is mandatory.

## [[StochSK] Block](@id stochsk_block)

!!! note

    The parameters in this block are valid for the `StochSK` solver only.

!!! warning

    If `solver = "StochSK"`, the `[StochSK]` block must be available.

### [method](@id stochsk_method)

*Definition:*

> How to determine the optimized ``\Theta`` parameter? The `StochSK` solver supports two different algorithms. They are
>
> * chi2min
> * chi2kink
>
> Usually, the `chi2min` algorithm is preferred. This algorithm is suggested by Shao and Sandvik *et al*. See [Stochastic Analytic Continuation 1](@ref san) for more details.

*Type:*

> String.

*Example:*

> method = "chi2min"

*Comment:*

> This parameter is mandatory.

### [nfine](@id stochsk_nfine)

*Definition:*

> Number of points of a very fine linear mesh. This mesh is for the ``\delta`` functions.

*Type:*

> Integer.

*Example:*

> nfine = 100000

*Comment:*

> This parameter is mandatory.

### [ngamm](@id stochsk_ngamm)

*Definition:*

> Number of ``\delta`` functions. Their superposition is used to mimic the spectral functions.

*Type:*

> Integer.

*Example:*

> ngamm = 1000

*Comment:*

> This parameter is mandatory.

### [nwarm](@id stochsk_nwarm)

*Definition:*

> Number of Monte Carlo thermalization steps.

*Type:*

> Integer.

*Example:*

> nwarm = 1000

*Comment:*

> This parameter is mandatory.

### [nstep](@id stochsk_nstep)

*Definition:*

> Number of Monte Carlo sweeping steps.

*Type:*

> Integer.

*Example:*

> nstep = 20000

*Comment:*

> This parameter is mandatory.

### [ndump](@id stochsk_ndump)

*Definition:*

> Intervals for monitoring Monte Carlo sweeps. For every `ndump` steps, the `StochSK` solver will try to output some useful information to help diagnosis.

*Type:*

> Integer.

*Example:*

> ndump = 200

*Comment:*

> This parameter is mandatory.

### [retry](@id stochsk_retry)

*Definition:*

> How often to recalculate the goodness-of-fit function (it is actually ``\chi^2``) to avoid numerical deterioration.

*Type:*

> Integer.

*Example:*

> retry = 10

*Comment:*

> This parameter is mandatory.

### [theta](@id stochsk_theta)

*Definition:*

> Starting value for the ``\Theta`` parameter. The `StochSK` solver always starts with a huge ``\Theta`` parameter, and then decreases it gradually.

*Type:*

> Float.

*Example:*

> theta = 1e+6

*Comment:*

> This parameter is mandatory.

### [ratio](@id stochsk_ratio)

*Definition:*

> Scaling factor for the ``\Theta`` parameter. It should be less than 1.0.

*Type:*

> Float.

*Example:*

> ratio = 0.9

*Comment:*

> This parameter is mandatory.

## [[StochOM] Block](@id stochom_block)

!!! note

    The parameters in this block are valid for the `StochOM` solver only.

!!! warning

    If `solver = "StochOM"`, the `[StochOM]` block must be available.

### [ntry](@id stochom_ntry)

*Definition:*

> Number of attempts to figure out the solution.

*Type:*

> Integer.

*Example:*

> ntry = 2000

*Comment:*

> This parameter is mandatory.

### [nstep](@id stochom_nstep)

*Definition:*

> Number of Monte Carlo steps per try.

*Type:*

> Integer.

*Example:*

> nstep = 1000

*Comment:*

> This parameter is mandatory.

### [nbox](@id stochom_nbox)

*Definition:*

> Number of boxes. Their superposition is used to construct the spectral functions.

*Type:*

> Integer.

*Example:*

> nbox = 100

*Comment:*

> This parameter is mandatory.

### [sbox](@id stochom_sbox)

*Definition:*

> Minimum area of the randomly generated boxes.

*Type:*

> Float.

*Example:*

> sbox = 0.005

*Comment:*

> This parameter is mandatory.

### [wbox](@id stochom_wbox)

*Definition:*

> Minimum width of the randomly generated boxes.

*Type:*

> Float.

*Example:*

> wbox = 0.02

*Comment:*

> This parameter is mandatory.

### [norm](@id stochom_norm)

*Definition:*

> Is the norm calculated? If `norm` is larger than 0.0, it denotes the normalization factor. If `norm` is smaller than 0.0, it means that the normalization condition is ignored.

*Type:*

> Float.

*Example:*

> norm = -1.0

*Comment:*

> This parameter is mandatory.

## [[StochPX] Block](@id stochpx_block)

!!! note

    The parameters in this block are valid for the `StochPX` solver only.

!!! warning

    If `solver = "StochPX"`, the `[StochPX]` block must be available.

!!! warning

    The `StochPX` solver is still in development. Please use it at your own risk.

### [method](@id stochpx_method)

*Definition:*

> How to evaluate the final spectral density? The `StochPX` solver supports two different algorithms. They are
>
> * mean
> * best
>
> If `method = "mean"`, then the solver will try to calculate an averaged spectrum from some selected `good` solutions. If `method = "best"`, then the solver will pick up the best solution (which should exhibit the smallest goodness-of-fit functional ``\chi^2``).

*Type:*

> String.

*Example:*

> method = "mean"

*Comment:*

> This parameter is mandatory. Note that the "mean" method is suitable for the condensed matter cases (broad and smooth peaks), while the "best" method is useful for the molecule cases (sharp peaks).

### [nfine](@id stochpx_nfine)

*Definition:*

> Number of grid points for a very fine mesh. This mesh is for the poles.

*Type:*

> Integer.

*Example:*

> nfine = 100000

*Comment:*

> This parameter is mandatory.

### [npole](@id stochpx_npole)

*Definition:*

> Number of poles on the real axis. These poles are used to mimic the Matsubara Green's function.

*Type:*

> Integer.

*Example:*

> npole = 200

*Comment:*

> This parameter is mandatory. For condensed matter cases, `npole` should be quite large. While for molecule cases, `npole` should be small.

### [ntry](@id stochpx_ntry)

*Definition:*

> Number of attempts to figure out the solution.

*Type:*

> Integer.

*Example:*

> ntry = 1000

*Comment:*

> This parameter is mandatory.

### [nstep](@id stochpx_nstep)

*Definition:*

> Number of Monte Carlo sweeping steps per attempt / try.

*Type:*

> Integer.

*Example:*

> nstep = 1000000

*Comment:*

> This parameter is mandatory. This parameter is related to the `npole` parameter. If `npole` is large, `nstep` could be small. If `npole` is small, `nstep` should be large.

### [theta](@id stochpx_theta)

*Definition:*

> Artificial inverse temperature ``\Theta``. When it is increased, the transition probabilities of Monte Carlo updates will decrease.

*Type:*

> Float.

*Example:*

> theta = 1e+6

*Comment:*

> This parameter is mandatory. The users can check the `stat.data` file to judge whether the `theta` parameter is reasonable.

### [eta](@id stochpx_eta)

*Definition:*

> Tiny distance from the real axis ``\eta``, which is used to reconstruct the retarded Green's function and the spectral density. When it is increased, the spectral density will be become more and more smooth.

*Type:*

> Float.

*Example:*

> eta = 1e-4

*Comment:*

> This parameter is mandatory.
