# Kernels

*Build kernel functions.*

The ACFlow toolkit supports six types of kernel functions. They are:

* Fermionic imaginary time kernel (`ktype = "fermi", grid = "ftime"`)
* Fermionic Matsubara kernel (`ktype = "fermi", grid = "ffreq"`)
* Bosonic imaginary time kernel (`ktype = "boson", grid = "btime"`)
* Bosonc Matsubara kernel (`ktype = "boson", grid = "bfreq"`)
* Symmetric bosonic imaginary time kernel (`ktype = "bsymm", grid = "btime"`)
* Symmetric bosonic Matsubara kernel (`ktype = "bsymm", grid = "bfreq"`)

Note that the `MaxEnt`, `StochAC`, and `StochSK` solvers rely on the `make_kernel()` function to provide the kernel function. However, the kernel function or matrix used in the `StochOM` and `StochPX` solvers are implemented in their own `calc_lambda()` functions.

## Contents

```@contents
Pages = ["kernel.md"]
Depth = 2
```

## Index

```@index
Pages = ["kernel.md"]
```

## Making Kernels

```@docs
build_kernel
build_kernel_symm
```

## Utilities

```@docs
make_blur
make_singular_space
make_gauss_peaks
```
