# Parameters

## Contents

```@contents
Pages = ["param.md"]
Depth = 2
```

## [BASE] block

!!! note

    The parameters in this block is valid for all the solvers.

### finput

*Definition:*

> Filename for the input data.

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
> Here, `MaxEnt` means the maximum entropy method, `StochAC` means the stochastic analytical continuation method (K. S. D. Beach's version), `StochSK` means the stochastic analytical continuation method (A. W. Sandvik's version), and `StochOM` means the stochastic optimization method.

*Type:*

> String.

*Example:*

> solver = "MaxEnt"

*Comment:*

> If `solver = "MaxEnt"`, then the `[MaxEnt]` block must be available in the configuration file. If `solver = "StochAC"`, then the `[StochAC]` block must be available in the configuration file. If `solver = "StochSK"`, then the `[StochSK]` block must be available in the configuration file. If `solver = "StochOM"`, then the `[StochOM]` block must be available in the configuration file.
>
> This parameter is mandatory.

### ktype

*Definition:*

> It denotes type of kernel function. Now **ACFlow.jl** supports three types of kernel functions. They are:
>
> * fermi
> * boson
> * bsymm
>
> Here, `fermi` means fermionic Kernel, `boson` means bosonic kernel, and `bsymm` means symmetric bosonic kernel. 

*Type:*

> String.

*Example:*

> ktype = "fermi"

*Comment:*

> This parameter is mandatory.

### mtype

*Definition:*

> Type of default model function. Now **ACFlow.jl** supports the following choices:
>
> * flat
> * gauss
> * 1gauss
> * 2gauss
> * lorentz
> * 1lorentz
> * 2lorentz
> * risedecay

*Type:*

> String.

*Example:*

> mtype = "flat"

*Comment:*

> This parameter is mandatory.

### grid

*Definition:* 

> This parameter specifies the grid's type for input data (imaginary axis). Now **ACFlow.jl** supports the following choices:
>
> * ftime
> * btime
> * ffreq
> * bfreq

*Type:*

> String.

*Example:*

> grid = "ftime"

*Comment:*

> This parameter is mandatory.

### mesh

*Definition:*

> This parameter specifies the mesh's type for output data (real axis). Now **ACFlow.jl** supports the following choices:
>
> * linear
> * tangent
> * lorentz
> * halflorentz

*Type:*

> String.

*Comment:*

> This parameter is mandatory.

### ngrid

*Definition:*

> Number of grid points.

*Type:*

> Integer.

*Example:*

> ngrid = 10

*Comment:*

> This parameter must be compatible with the input data.

### nmesh

*Definition:*

> Number of mesh points.

*Type:*

> Integer.

*Example:*

> nmesh = 501

*Comment:*

> This parameter is mandatory.

### wmax

*Definition:*

> Right boundary (maximum value) of mesh.

*Type:*

> Float.

*Example:*

> wmax = 10.0

*Comment:*

> This parameter is mandatory.

### wmin

*Definition:*

> Left boundary (minimum value) of mesh.

*Type:*

> Float.

*Example:*

> wmin = -10.0

*Comment:*

> This parameter is mandatory.

### beta

*Definition:*

> Inverse temperature . It is equal to 1/T.

*Type:*

> Float.

*Example:*

> beta = 10.0

*Comment:*

> This parameter must be compatible with the input data and grid.

### offdiag

*Definition:*

*Type:*

*Comment:*

### pmodel

*Definition:*

*Type:*

*Comment:*

### pmesh

*Definition:*

*Type:*

*Comment:*

### exclude

*Definition:*

*Type:*

*Comment:*

## [MaxEnt] block

### method

*Definition:*

*Type:*

*Comment:*

### nalph

*Definition:* Total number of the chosen α parameters.

*Type:* Integer.

*Comment:*

### alpha

*Definition:* Starting value for the α parameter.

*Type:* Float.

*Comment:* It should be a very large number, such as 1e9 - 1e13.

### ratio

*Definition:*

*Type:*

*Comment:*

### blur

*Definition:*

*Type:*

*Comment:*

## [StochAC] block

### nfine

*Definition:*

*Type:*

*Comment:*

### ngamm

*Definition:* Number of δ functions, which is used to mimic the spectral functions.

*Type:* Integer.

*Comment:*

### nwarm

*Definition:*

*Type:*

*Comment:*

### nstep

*Definition:*

*Type:*

*Comment:*

### ndump

*Definition:*

*Type:*

*Comment:*

### nalph

*Definition:*

*Type:*

*Comment:*

### alpha

*Definition:*

*Type:*

*Comment:*

### ratio

*Definition:*

*Type:*

*Comment:*

## [StochSK] block

### nfine

*Definition:*

*Type:*

*Comment:*

### ngamm

*Definition:*

*Type:*

*Comment:*

### nwarm

*Definition:*

*Type:*

*Comment:*

### nstep

*Definition:*

*Type:*

*Comment:*

### ndump

*Definition:*

*Type:*

*Comment:*

### retry

*Definition:*

*Type:*

*Comment:*

### theta

*Definition:*

*Type:*

*Comment:*

### ratio

*Definition:*

*Type:*

*Comment:*

## [StochOM] block

### ntry

*Definition:*

*Type:*

*Comment:*

### nstep

*Definition:*

*Type:*

*Comment:*

### nbox

*Definition:*

*Type:*

*Comment:*

### sbox

*Definition:*

*Type:*

*Comment:*

### wbox

*Definition:*

*Type:*

*Comment:*

### norm

*Definition:*

*Type:*

*Comment:*
