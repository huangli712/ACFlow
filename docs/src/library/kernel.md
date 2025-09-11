*Build kernel functions.*

The ACFlow toolkit supports twelve types of kernel functions. They are:

* Fermionic imaginary time kernel (`ktype = "fermi", grid = "ftime"`)
* Fermionic fragment imaginary time kernel (`ktype = "fermi", grid = "fpart"`)
* Fermionic Matsubara kernel (`ktype = "fermi", grid = "ffreq"`)
* Fermionic fragment Matsubara kernel (`ktype = "fermi", grid = "ffrag"`)
* Bosonic imaginary time kernel (`ktype = "boson", grid = "btime"`)
* Bosonic fragment imaginary time kernel (`ktype = "boson", grid = "bpart"`)
* Bosonc Matsubara kernel (`ktype = "boson", grid = "bfreq"`)
* Bosonc fragment Matsubara kernel (`ktype = "boson", grid = "bfrag"`)
* Symmetric bosonic imaginary time kernel (`ktype = "bsymm", grid = "btime"`)
* Symmetric bosonic fragment imaginary time kernel (`ktype = "bsymm", grid = "bpart"`)
* Symmetric bosonic Matsubara kernel (`ktype = "bsymm", grid = "bfreq"`)
* Symmetric bosonic fragment Matsubara kernel (`ktype = "bsymm", grid = "bfrag"`)

Note that the `MaxEnt`, `StochAC`, and `StochSK` solvers rely on the `make_kernel()` function to provide the kernel function. However, the kernel function or matrix used in the `StochOM` and `StochPX` solvers are implemented in their own `calc_lambda()` functions.

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
