# Parameters

## Contents

```@contents
Pages = ["param.md"]
Depth = 3
```

## [BASE] block

!!! note

    The parameters in this block is valid for all the solvers.

### finput

*Definition:*

> Filename for the input data. The data should be saved in a column-wised and formated (CSV-like) text file.

*Type:*

> String.

*Example:*

> finput = "green.data"

*Comment:*

> This parameter is mandatory.

### solver

*Definition:*

> This parameter specifies the solvers that used to solve the analytical continuation problem. Now **ACFlow.jl** supports four different solvers. They are as follows:
>
> * MaxEnt
> * StochAC
> * StochSK
> * StochOM
>
> Here, `MaxEnt` means the maximum entropy method. The `MaxEnt` solver can be used to treat the correlators in Matsubara frequency or imaginary-time axis. If `solver = "MaxEnt"`, then the `[MaxEnt]` block must be available in the configuration file.
>
>`StochAC` means the stochastic analytical continuation method (K. S. D. Beach's version). The `StochAC` solver can be used to treat the correlators in imaginary-time axis only. If `solver = "StochAC"`, then the `[StochAC]` block must be available in the configuration file.
>
> `StochSK` means the stochastic analytical continuation method (A. W. Sandvik's version). The `StochSK` solver can be used to treat the correlators in imaginary-time axis only. If `solver = "StochSK"`, then the `[StochSK]` block must be available in the configuration file.
>
> `StochOM` means the stochastic optimization method. The `StochOM` solver can be used to treat the correlators in Matsubara frequency axis only. If `solver = "StochOM"`, then the `[StochOM]` block must be available in the configuration file.

*Type:*

> String.

*Example:*

> solver = "MaxEnt"

*Comment:*

> This parameter is mandatory.

### ktype

*Definition:*

> It denotes the type of kernel functions. Now **ACFlow.jl** supports three types of kernel functions. They are:
>
> * fermi
> * boson
> * bsymm
>
> Here, `fermi` means fermionic kernel, `boson` means bosonic kernel, and `bsymm` means symmetric bosonic kernel. 

*Type:*

> String.

*Example:*

> ktype = "fermi"

*Comment:*

> This parameter is mandatory.

### mtype

*Definition:*

> It denotes the type of default model functions. Now **ACFlow.jl** supports the following choices:
>
> * flat
> * gauss
> * 1gauss
> * 2gauss
> * lorentz
> * 1lorentz
> * 2lorentz
> * risedecay
>
> Here, `flat` means the flat model (i.e, constant), `gauss` means the Gaussian model, `1gauss` means the Shifted Gaussian model, `2gauss` means the Two Gaussian model, `lorentz` means the Lorentzian model, `1lorentz` means the Shifted Lorentzian model, `2lorentz` means the Two Lorentzian model, and `risedecay` means the Rise-And-Decay model. As for detailed formula for these models, please refer to `src/model.jl`.

*Type:*

> String.

*Example:*

> mtype = "flat"

*Comment:*

> This parameter is mandatory. Only the `MaxEnt` solver need these model functions. The `StochAC` solver only supports the `flat` model. The `StochSK` and `StochOM` are free of model functions.

### grid

*Definition:* 

> This parameter specifies the grid's type for input data in imaginary axis. Now **ACFlow.jl** supports the following choices:
>
> * ftime
> * btime
> * ffreq
> * bfreq
>
> Here, `ftime` means fermionic and imaginary-time, `btime` means bosonic and imaginary-time, `ffreq` means fermionic and Matsubara frequency, and `bfreq` means bosonic and Matsubara frequency.

*Type:*

> String.

*Example:*

> grid = "ftime"

*Comment:*

> This parameter is mandatory.

### mesh

*Definition:*

> This parameter specifies the mesh's type for output data (usually the spectral functions) in real axis. Now **ACFlow.jl** supports the following choices:
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

> This parameter is mandatory.

### ngrid

*Definition:*

> Number of grid points. The parameter, together with the `beta` and `grid` parameters, control the generation of grid for input data. 

*Type:*

> Integer.

*Example:*

> ngrid = 10

*Comment:*

> This parameter must be compatible with the input data.

### nmesh

*Definition:*

> Number of mesh points. The parameter, together with the `wmax`, `wmin`, and `mesh` parameters, control the generation of mesh for output data. 

*Type:*

> Integer.

*Example:*

> nmesh = 501

*Comment:*

> This parameter is mandatory.

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

### beta

*Definition:*

> Inverse temperature ``\beta``. It is equal to 1/T.

*Type:*

> Float.

*Example:*

> beta = 10.0

*Comment:*

> This parameter must be compatible with the input data and grid.

### offdiag

*Definition:*

> Is it the offdiagonal part in matrix-valued function.

*Type:*

> Bool.

*Example:*

> offdiag = false

*Comment:*

> This parameter is useful for the `MaxEnt` solver only.

### pmodel

*Definition:*

> Additional parameters for customizing the model.

*Type:*

> Array.

*Example:*

>

*Comment:*

> This parameter is optional.

### pmesh

*Definition:*

> Additional parameters for customizing the mesh.

*Type:*

> Array.

*Example:*

>

*Comment:*

> This parameter is optional.

### exclude

*Definition:*

> Restriction of the energy range of the spectrum.

*Type:*

> Array.

*Example:*

> exclude = [[-8.0,-4.0],[4.0,8.0]]

*Comment:*

> This parameter is optional.

## [MaxEnt] block

### method

*Definition:*

> How to determine the optimized α parameter? The `MaxEnt` solver supports four different algorithms. They are
>
> * historic
> * classic
> * bryan
> * chi2kink

*Type:*

> String.

*Example:*

> method = "bryan"

*Comment:*

> This parameter is mandatory.

### nalph

*Definition:*

> Total number of the chosen α parameters.

*Type:*

> Integer.

*Example:*

> nalph = 12

*Comment:*

> This parameter is mandatory.

### alpha

*Definition:*

> Starting value for the α parameter.

*Type:*

> Float.

*Example:*

> alpha = 1e9

*Comment:*

> It should be a very large number, such as 1e9 - 1e13.

### ratio

*Definition:*

> Scaling factor for the α parameter.

*Type:*

> Float.

*Example:*

>ratio = 10.0

*Comment:*

> This parameter is mandatory.

### blur

*Definition:*

> Shall we preblur the kernel and spectrum?

*Type:*

> Float.

*Example:*

> blur = -1.0

*Comment:*

> This parameter is mandatory.

## [StochAC] block

### nfine

*Definition:*

> Number of points of a very fine linear mesh.

*Type:*

> Integer.

*Example:*

> nfine = 10000

*Comment:*

> This parameter is mandatory.

### ngamm

*Definition:* 

> Number of δ functions, which is used to mimic the spectral functions.

*Type:*

> Integer.

*Example:*

> ngamm = 512

*Comment:*

> This parameter is mandatory.

### nwarm

*Definition:*

> Number of Monte Carlo warmup steps.

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

> Intervals for monitoring Monte Carlo sweeps.

*Type:*

> Integer.

*Example:*

> ndump = 40000

*Comment:*

> This parameter is mandatory.

### nalph

*Definition:*

> Total number of the chosen α parameters.

*Type:*

> Integer.

*Example:*

> nalph = 20

*Comment:*

> This parameter is mandatory.

### alpha

*Definition:*

> Starting value for the α parameter.

*Type:*

> Float.

*Example:*

> alpha = 1.0

*Comment:*

> This parameter is mandatory.

### ratio

*Definition:*

> Scaling factor for the α parameter.

*Type:*

> Float.

*Example:*

> ratio = 1.2

*Comment:*

> This parameter is mandatory.

## [StochSK] block

### nfine

*Definition:*

> Number of points of a very fine linear mesh.

*Type:*

> Integer.

*Example:*

> nfine = 100000

*Comment:*

> This parameter is mandatory.

### ngamm

*Definition:*

> Number of δ functions.

*Type:*

> Integer.

*Example:*

> ngamm = 1000

*Comment:*

> This parameter is mandatory.

### nwarm

*Definition:*

> Number of Monte Carlo warmup steps.

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

> Intervals for monitoring Monte Carlo sweeps.

*Type:*

> Integer.

*Example:*

> ndump = 200

*Comment:*

> This parameter is mandatory.

### retry

*Definition:*

> How often to recalculate the goodness function.

*Type:*

> Integer.

*Example:*

> retry = 10

*Comment:*

> This parameter is mandatory.

### theta

*Definition:*

> Starting value for the Θ parameter.

*Type:*

> Float.

*Example:*

> theta = 1e+6

*Comment:*

> This parameter is mandatory.

### ratio

*Definition:*

> Scaling factor for the Θ parameter.

*Type:*

> Float.

*Example:*

> ratio = 0.9

*Comment:*

> This parameter is mandatory.

## [StochOM] block

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

> Number of boxes to construct the spectrum.

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

> Is the norm calculated.

*Type:*

> Float.

*Example:*

> norm = -1.0

*Comment:*

> This parameter is mandatory.
