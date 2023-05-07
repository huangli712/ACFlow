# [Parameters](@id param)

*A comprehensive dictionary about parameters*

The official configuration file for the ACFlow toolkit is `case.toml`. This page contains all the valid parameters that can appear in `case.toml`. As for the format of `case.toml`, please look at [`case.toml`](input.md).

## Contents

```@contents
Pages = ["param.md"]
Depth = 3
```

## [BASE] Block

!!! note

    This block is mandatory. The parameters in this block is useful for all the solvers.

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

> This parameter specifies the solvers that used to solve the analytical continuation problem. Now the ACFlow toolkit supports five different solvers. They are as follows:
>
> * MaxEnt
> * StochAC
> * StochSK
> * StochOM
> * StochPX
>
> Here, `MaxEnt` means the maximum entropy method. The `MaxEnt` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "MaxEnt"`, then the `[MaxEnt]` block must be available in the configuration file.
>
>`StochAC` means the stochastic analytical continuation method (K. S. D. Beach's algorithm). The `StochAC` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochAC"`, then the `[StochAC]` block must be available in the configuration file.
>
> `StochSK` means the stochastic analytical continuation method (A. W. Sandvik's algorithm). The `StochSK` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochSK"`, then the `[StochSK]` block must be available in the configuration file.
>
> `StochOM` means the stochastic optimization method. The `StochOM` solver can be used to treat the correlators in Matsubara frequency or imaginary time axis. If `solver = "StochOM"`, then the `[StochOM]` block must be available in the configuration file.
>
> `StochPX` means the stochastic pole expansion method. The `StochPX` solver can be used to treat the correlators in Matsubara frequency axis only. If `solver = "StochPX"`, then the `[StochPX]` block must be available in the configuration file.

!!! warning

    For the `StochOM` solver, if the correlators are defined in imaginary time axis, they must be bosonic. In other words, the `StochOM` solver does not support analytical continuation of fermionic imaginary time correlation function.

*Type:*

> String.

*Example:*

> solver = "MaxEnt"

*Comment:*

> This parameter is mandatory.

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
> Especially, if `mtype = "file"`, then the default model function is encoded in `model.inp`. ACFlow will read this file and initialize the default model function automatically. Be careful, the mesh for this model function must be consistent with the one used in the analytical continuation calculations.
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
> * btime
> * fpart
> * bpart
> * ffreq
> * bfreq
>
> Here, `ftime` means fermionic and imaginary time, `btime` means bosonic and imaginary time, `ffreq` means fermionic and Matsubara frequency, and `bfreq` means bosonic and Matsubara frequency. `fpart` means fermionic and imaginary time as well, but the grid of imaginary time might be incomplete. `bpart` is similar to `fpart`, but it is for the bosonic case.

*Type:*

> String.

*Example:*

> grid = "ftime"

*Comment:*

> This parameter is mandatory. It must be compatible with the `ktype` parameter. See also [`ngrid`](@ref ngrid).

!!! warning

    If the `StochOM` solver is employed, the `grid` parameter should not be "ftime" or "fpart".

!!! warning

    If the `StochPX` solver is employed, the `grid` parameter should be "ffreq" or "bfreq".

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

### [nmesh](@id nmesh)

*Definition:*

> Number of mesh points. The parameter, together with the `wmax`, `wmin`, and `mesh` parameters, controls the generation of mesh for output data.

*Type:*

> Integer.

*Example:*

> nmesh = 501

*Comment:*

> This parameter is mandatory. See also [`mesh`](@ref mesh).

### wmax

*Definition:*

> Right boundary (maximum value) of mesh. Note that `wmax` should be always greater than `wmin`.

*Type:*

> Float.

*Example:*

> wmax = 10.0

*Comment:*

> This parameter is mandatory.

### wmin

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

### beta

*Definition:*

> Inverse temperature ``\beta``. It is equal to ``1/T``.

*Type:*

> Float.

*Example:*

> beta = 10.0

*Comment:*

> This parameter is mandatory. This parameter must be compatible with the input data and grid. Specifically, for the imaginary time axis, the last grid point should be ``\beta``. As for the Matsubara frequency axis, the difference between two successive grid points should be ``\pi/\beta``.

### offdiag

*Definition:*

> Is the input correlator the offdiagonal part in matrix-valued function? As for the offdiagonal correlator, the corresponding spectral function might be not positive-definite. Some tricks have been implemented to cure this issue.

*Type:*

> Bool.

*Example:*

> offdiag = false

*Comment:*

> This parameter is mandatory. This parameter is useful for the `MaxEnt` solver only.

!!! warning

    Now only the `MaxEnt` solver supports this parameter.

### fwrite

*Definition:*

> Are the analytical continuation results written into external files? If it is false, then only the terminal output is retained and all the other outputs are disable. By default (if this parameter is missing or true), the files should be generated.

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

### exclude

*Definition:*

> Restriction of the energy range of the calculated spectral functions. This features is implemented by the `StochAC`, `StochSK`, `StochOM`, and `StochPX` solvers. In these solvers, the ``\delta`` or `box` functions, which are used to mimic the spectral functions, are restricted to live out of the given energy ranges. For example, `exclude = [8.0,16.0]` means that the energy range `[8.0,16.0]` is strictly forbidden.

*Type:*

> Array.

*Example:*

> exclude = [[-8.0,-4.0],[4.0,8.0]]

*Comment:*

> This parameter is optional. If you are using the `MaxEnt` solver, this parameter will be ignored.

## [MaxEnt] Block

!!! note

    The parameters in this block is valid for the `MaxEnt` solver only.

!!! warning

    If `solver = "Maxent"`, the `[MaxEnt]` block must be available.

### method

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

### stype

*Definition:*

> Type of the entropic factor. The `MaxEnt` solver supports two schemes. They are
>
> * sj
> * br
>
> Here, `sj` means the Shannon-Jaynes entropy, while `br` means the Bayesian Reconstruction entropy. Usually, the Shannon-Jaynes entropy is preferred, since with it the positivity of the generated spectrum is always guaranteed.

*Type:*

> String.

*Example:*

> stype = "sj"

*Comment:*

> This parameter is mandatory. As for the underlying principles of these entropic factors, please see [Maximum Entropy Method](@ref mem).

### nalph

*Definition:*

> Total number of the chosen ``\alpha`` parameters.

*Type:*

> Integer.

*Example:*

> nalph = 12

*Comment:*

> This parameter is mandatory. Only the `chi2kink` algorithm needs this parameter to control the number of ``\alpha`` parameters.

### alpha

*Definition:*

> Starting value for the ``\alpha`` parameter. The `MaxEnt` solver always starts with a huge ``\alpha`` parameter, and then decreases it gradually.

*Type:*

> Float.

*Example:*

> alpha = 1e9

*Comment:*

> This parameter is mandatory. It should be a very large number, such as ``10^9 \sim 10^{13}``.

### ratio

*Definition:*

> Scaling factor for the ``\alpha`` parameter. The next ``\alpha`` is equal to the current ``\alpha`` divided by `ratio`.

*Type:*

> Float.

*Example:*

>ratio = 10.0

*Comment:*

> This parameter is mandatory. It muse be larger than 1.0.

### blur

*Definition:*

> Sometimes, the kernel functions and spectral functions can be preblurred to obtain smoother results. Shall we preblur them? If `blur` is larger than zero, then it means the blur parameter. If `blur` is smaller than zero, then it means that the preblur feature is disable.

*Type:*

> Float.

*Example:*

> blur = -1.0

*Comment:*

> This parameter is mandatory.

## [StochAC] Block

!!! note

    The parameters in this block is valid for the `StochAC` solver only.

!!! warning

    If `solver = "StochAC"`, the `[StochAC]` block must be available.

### nfine

*Definition:*

> Number of points of a very fine linear mesh. This mesh is for the ``\delta`` functions.

*Type:*

> Integer.

*Example:*

> nfine = 10000

*Comment:*

> This parameter is mandatory.

### ngamm

*Definition:*

> Number of ``\delta`` functions. Their superposition is used to mimic the spectral functions.

*Type:*

> Integer.

*Example:*

> ngamm = 512

*Comment:*

> This parameter is mandatory.

### nwarm

*Definition:*

> Number of Monte Carlo thermalization steps.

*Type:*

> Integer.

*Example:*

> nwarm = 4000

*Comment:*

> This parameter is mandatory.

### nstep

*Definition:*

> Number of Monte Carlo sweeping steps.

*Type:*

> Integer.

*Example:*

> nstep = 4000000

*Comment:*

> This parameter is mandatory.

### ndump

*Definition:*

> Intervals for monitoring Monte Carlo sweeps. For every `ndump` steps, the `StochAC` solver will try to output some useful information to help diagnosis.

*Type:*

> Integer.

*Example:*

> ndump = 40000

*Comment:*

> This parameter is mandatory.

### nalph

*Definition:*

> Total number of the chosen ``\alpha`` parameters.

*Type:*

> Integer.

*Example:*

> nalph = 20

*Comment:*

> This parameter is mandatory.

### alpha

*Definition:*

> Starting value for the ``\alpha`` parameter. The `StochAC` solver always starts with a small ``\alpha`` parameter, and then increases it gradually.

*Type:*

> Float.

*Example:*

> alpha = 1.0

*Comment:*

> This parameter is mandatory.

### ratio

*Definition:*

> Scaling factor for the ``\alpha`` parameter. It should be larger than 1.0.

*Type:*

> Float.

*Example:*

> ratio = 1.2

*Comment:*

> This parameter is mandatory.

## [StochSK] Block

!!! note

    The parameters in this block is valid for the `StochSK` solver only.

!!! warning

    If `solver = "StochSK"`, the `[StochSK]` block must be available.

### method

*Definition:*

> How to determine the optimized ``\Theta`` parameter? The `StochSK` solver supports two different algorithms. They are
>
> * chi2min
> * chi2kink
>
> Usually, the `chi2min` algorithm is preferred. This algorithm is suggested by Shao and Sandvik *et al*. See [Stochastic Analytical Continuation 1](@ref san) for more details.

*Type:*

> String.

*Example:*

> method = "chi2min"

*Comment:*

> This parameter is mandatory.

### nfine

*Definition:*

> Number of points of a very fine linear mesh. This mesh is for the ``\delta`` functions.

*Type:*

> Integer.

*Example:*

> nfine = 100000

*Comment:*

> This parameter is mandatory.

### ngamm

*Definition:*

> Number of ``\delta`` functions. Their superposition is used to mimic the spectral functions.

*Type:*

> Integer.

*Example:*

> ngamm = 1000

*Comment:*

> This parameter is mandatory.

### nwarm

*Definition:*

> Number of Monte Carlo thermalization steps.

*Type:*

> Integer.

*Example:*

> nwarm = 1000

*Comment:*

> This parameter is mandatory.

### nstep

*Definition:*

> Number of Monte Carlo sweeping steps.

*Type:*

> Integer.

*Example:*

> nstep = 20000

*Comment:*

> This parameter is mandatory.

### ndump

*Definition:*

> Intervals for monitoring Monte Carlo sweeps. For every `ndump` steps, the `StochSK` solver will try to output some useful information to help diagnosis.

*Type:*

> Integer.

*Example:*

> ndump = 200

*Comment:*

> This parameter is mandatory.

### retry

*Definition:*

> How often to recalculate the goodness-of-fit function (it is actually ``\chi^2``) to avoid numerical deterioration.

*Type:*

> Integer.

*Example:*

> retry = 10

*Comment:*

> This parameter is mandatory.

### theta

*Definition:*

> Starting value for the ``\Theta`` parameter. The `StochSK` solver always starts with a huge ``\Theta`` parameter, and then decreases it gradually.

*Type:*

> Float.

*Example:*

> theta = 1e+6

*Comment:*

> This parameter is mandatory.

### ratio

*Definition:*

> Scaling factor for the ``\Theta`` parameter. It should be less than 1.0.

*Type:*

> Float.

*Example:*

> ratio = 0.9

*Comment:*

> This parameter is mandatory.

## [StochOM] Block

!!! note

    The parameters in this block is valid for the `StochOM` solver only.

!!! warning

    If `solver = "StochOM"`, the `[StochOM]` block must be available.

### ntry

*Definition:*

> Number of attempts to figure out the solution.

*Type:*

> Integer.

*Example:*

> ntry = 2000

*Comment:*

> This parameter is mandatory.

### nstep

*Definition:*

> Number of Monte Carlo steps per try.

*Type:*

> Integer.

*Example:*

> nstep = 1000

*Comment:*

> This parameter is mandatory.

### nbox

*Definition:*

> Number of boxes. Their superposition is used to construct the spectral functions.

*Type:*

> Integer.

*Example:*

> nbox = 100

*Comment:*

> This parameter is mandatory.

### sbox

*Definition:*

> Minimum area of the randomly generated boxes.

*Type:*

> Float.

*Example:*

> sbox = 0.005

*Comment:*

> This parameter is mandatory.

### wbox

*Definition:*

> Minimum width of the randomly generated boxes.

*Type:*

> Float.

*Example:*

> wbox = 0.02

*Comment:*

> This parameter is mandatory.

### norm

*Definition:*

> Is the norm calculated? If `norm` is larger than 0.0, it denotes the normalization factor. If `norm` is smaller than 0.0, it means that the normalization condition is ignored.

*Type:*

> Float.

*Example:*

> norm = -1.0

*Comment:*

> This parameter is mandatory.

## [StochPX] Block

!!! note

    The parameters in this block is valid for the `StochPX` solver only.

!!! warning

    If `solver = "StochPX"`, the `[StochPX]` block must be available.

!!! warning

    The `StochPX` solver is still in development. Please use it at your own risk.

### method

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

### nfine

*Definition:*

> Number of points of a very fine linear mesh. This mesh is for the poles.

*Type:*

> Integer.

*Example:*

> nfine = 100000

*Comment:*

> This parameter is mandatory.

### npole

*Definition:*

> Number of poles on the real axis, which is used to mimic the Matsubara Green's function.

*Type:*

> Integer.

*Example:*

> npole = 200

*Comment:*

> This parameter is mandatory. For condensed matter cases, `npole` should be quite large. While for molecule cases, `npole` should be small.

### ntry

*Definition:*

> Number of attempts to figure out the solution.

*Type:*

> Integer.

*Example:*

> ntry = 1000

*Comment:*

> This parameter is mandatory.

### nstep

*Definition:*

> Number of Monte Carlo sweeping steps per attempt / try.

*Type:*

> Integer.

*Example:*

> nstep = 1000000

*Comment:*

> This parameter is mandatory. This parameter is related to the `npole` parameter. If `npole` is large, `nstep` could be small. If `npole` is small, `nstep` should be large.

### theta

*Definition:*

> Artificial inverse temperature ``\Theta``. When it is increased, the transition probabilities of Monte Carlo updates will decrease.

*Type:*

> Float.

*Example:*

> theta = 1e+6

*Comment:*

> This parameter is mandatory. The users can check the `stat.data` file to judge whether the `theta` parameter is reasonable.

### eta

*Definition:*

> Tiny distance from the real axis ``\eta``, which is used to reconstruct the retarded Green's function and the spectral density. When it is increased, the spectral density will be become more and more smooth.

*Type:*

> Float.

*Example:*

> eta = 1e-4

*Comment:*

> This parameter is mandatory.
