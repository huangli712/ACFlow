#
# Project : Gardenia
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/01/24
#

#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

#=
### *Customized Dictionaries*
=#

#=
*Remarks* :

The values in the following dictionaries are actually arrays, which
usually contain four elements:
* Element[1] -> Actually value.
* Element[2] -> If it is 1, this key-value pair is mandatory.
                If it is 0, this key-value pair is optional.
* Element[3] -> Numerical type (A julia Symbol).
* Element[4] -> Brief explanations.

The following dictionaries are used as global variables.
=#

"""
    PBASE

Dictionary for configuration parameters: general setup.
"""
const PBASE    = Dict{String,ADT}(
    "finput"  => [missing, 1, :String, "Filename for input data"],
    "solver"  => [missing, 1, :String, "Solver for the analytic continuation problem"],
    "ktype"   => [missing, 1, :String, "Type of kernel function"],
    "mtype"   => [missing, 1, :String, "Type of default model function"],
    "grid"    => [missing, 1, :String, "Grid for input data (imaginary axis)"],
    "mesh"    => [missing, 1, :String, "Mesh for output data (real axis)"],
    "ngrid"   => [missing, 1, :I64   , "Number of grid points"],
    "nmesh"   => [missing, 1, :I64   , "Number of mesh points"],
    "wmax"    => [missing, 1, :F64   , "Right boundary (maximum value) of output mesh"],
    "wmin"    => [missing, 1, :F64   , "Left boundary (minimum value) of output mesh"],
    "beta"    => [missing, 1, :F64   , "Inverse temperature"],
    "offdiag" => [missing, 1, :Bool  , "Is it the offdiagonal part in matrix-valued function"],
    "fwrite"  => [missing, 0, :Bool  , "Are the analytic continuation results written into files"],
    "pmodel"  => [missing, 0, :Array , "Additional parameters for customizing the model"],
    "pmesh"   => [missing, 0, :Array , "Additional parameters for customizing the mesh"],
    "exclude" => [missing, 0, :Array , "Restriction of the energy range of the spectrum"],
)

# Default parameters for PBASE
const _PBASE   = Dict{String,Any}(
    "finput"  => "green.data",
    "solver"  => "MaxEnt",
    "ktype"   => "fermi",
    "mtype"   => "flat",
    "grid"    => "ffreq",
    "mesh"    => "linear",
    "ngrid"   => 10,
    "nmesh"   => 501,
    "wmax"    => 5.0,
    "wmin"    => -5.0,
    "beta"    => 10.0,
    "offdiag" => false,
)

"""
    PMaxEnt

Dictionary for configuration parameters:
the maximum entropy method.
"""
const PMaxEnt  = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized α parameter"],
    "stype"   => [missing, 1, :String, "Type of the entropy term"],
    "nalph"   => [missing, 1, :I64   , "Total number of the chosen α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we preblur the kernel and spectrum"],
)

# Default parameters for PMaxEnt
const _PMaxEnt = Dict{String,Any}(
    "method"  => "chi2kink",
    "stype"   => "sj",
    "nalph"   => 12,
    "alpha"   => 1e9,
    "ratio"   => 10.0,
    "blur"    => -1.0, # Negative value means off.
)

"""
    PBarRat

Dictionary for configuration parameters:
the barycentric rational function approximation method.
"""
const PBarRat  = Dict{String,ADT}(
    "atype"   => [missing, 1, :String, "Possible type of the spectrum"],
    "denoise" => [missing, 1, :String, "How to denoise the input data"],
    "epsilon" => [missing, 1, :F64   , "Threshold for the Prony approximation"],
    "pcut"    => [missing, 1, :F64   , "Cutoff for unphysical poles"],
    "eta"     => [missing, 1, :F64   , "Tiny distance from the real axis"],
)

# Default parameters for PBarRat
const _PBarRat = Dict{String,Any}(
    "atype"   => "cont",
    "denoise" => "prony",
    "epsilon" => 1e-10,
    "pcut"    => 1e-3,
    "eta"     => 1e-2,
)

"""
    PNevanAC

Dictionary for configuration parameters:
the Nevanlinna analytical continuation method.
"""
const PNevanAC = Dict{String,ADT}(
    "pick"    => [missing, 1, :Bool  , "Check the Pick criterion or not"],
    "hardy"   => [missing, 1, :Bool  , "Perform Hardy basis optimization or not"],
    "hmax"    => [missing, 1, :I64   , "Upper cut off of Hardy order"],
    "alpha"   => [missing, 1, :F64   , "Regulation parameter for smooth norm"],
    "eta"     => [missing, 1, :F64   , "Tiny distance from the real axis"],
)

# Default parameters for PNevanAC
const _PNevanAC= Dict{String,Any}(
    "pick"    => true,
    "hardy"   => true,
    "hmax"    => 50,
    "alpha"   => 1e-4,
    "eta"     => 1e-2,
)

"""
    PStochAC

Dictionary for configuration parameters:
the stochastic analytic continuation method (K. S. D. Beach's version).
"""
const PStochAC = Dict{String,ADT}(
    "nfine"   => [missing, 1, :I64   , "Number of points of a very fine linear mesh"],
    "ngamm"   => [missing, 1, :I64   , "Number of δ functions"],
    "nwarm"   => [missing, 1, :I64   , "Number of Monte Carlo thermalization steps"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo sweeping steps"],
    "ndump"   => [missing, 1, :I64   , "Intervals for monitoring Monte Carlo sweeps"],
    "nalph"   => [missing, 1, :I64   , "Total number of the chosen α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
)

# Default parameters for PStochAC
const _PStochAC= Dict{String,Any}(
    "nfine"   => 10000,
    "ngamm"   => 512,
    "nwarm"   => 4000,
    "nstep"   => 4000000,
    "ndump"   => 40000,
    "nalph"   => 20,
    "alpha"   => 1.00,
    "ratio"   => 1.20,
)

"""
    PStochSK

Dictionary for configuration parameters:
the stochastic analytic continuation method (A. W. Sandvik's version).
"""
const PStochSK = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized Θ parameter"],
    "nfine"   => [missing, 1, :I64   , "Number of points of a very fine linear mesh"],
    "ngamm"   => [missing, 1, :I64   , "Number of δ functions"],
    "nwarm"   => [missing, 1, :I64   , "Number of Monte Carlo thermalization steps"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo sweeping steps"],
    "ndump"   => [missing, 1, :I64   , "Intervals for monitoring Monte Carlo sweeps"],
    "retry"   => [missing, 1, :I64   , "How often to recalculate the goodness function"],
    "theta"   => [missing, 1, :F64   , "Starting value for the Θ parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the Θ parameter"],
)

# Default parameters for PStochSK
const _PStochSK= Dict{String,Any}(
    "method"  => "chi2min",
    "nfine"   => 100000,
    "ngamm"   => 1000,
    "nwarm"   => 1000,
    "nstep"   => 20000,
    "ndump"   => 200,
    "retry"   => 10,
    "theta"   => 1e+6,
    "ratio"   => 0.90,
)

"""
    PStochOM

Dictionary for configuration parameters:
the stochastic optimization method.
"""
const PStochOM = Dict{String,ADT}(
    "ntry"    => [missing, 1, :I64   , "Number of attempts (tries) to seek the solution"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo steps per attempt / try"],
    "nbox"    => [missing, 1, :I64   , "Number of boxes to construct the spectrum"],
    "sbox"    => [missing, 1, :F64   , "Minimum area of the randomly generated boxes"],
    "wbox"    => [missing, 1, :F64   , "Minimum width of the randomly generated boxes"],
    "norm"    => [missing, 1, :F64   , "Is the norm calculated"],
)

# Default parameters for PStochOM
const _PStochOM= Dict{String,Any}(
    "ntry"    => 2000,
    "nstep"   => 1000,
    "nbox"    => 100,
    "sbox"    => 0.005,
    "wbox"    => 0.02,
    "norm"    => -1.0, # Negative value means off.
)

"""
    PStochPX

Dictionary for configuration parameters:
the stochastic pole expansion method.
"""
const PStochPX = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to evaluate the final spectral density"],
    "nfine"   => [missing, 1, :I64   , "Number of grid points for a very fine mesh"],
    "npole"   => [missing, 1, :I64   , "Number of poles"],
    "ntry"    => [missing, 1, :I64   , "Number of attempts (tries) to seek the solution"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo steps per attempt / try"],
    "theta"   => [missing, 1, :F64   , "Artificial inverse temperature"],
    "eta"     => [missing, 1, :F64   , "Tiny distance from the real axis"],
)

# Default parameters for PStochPX
const _PStochPX= Dict{String,Any}(
    "method"  => "mean",
    "nfine"   => 100000,
    "npole"   => 200,
    "ntry"    => 1000,
    "nstep"   => 1000000,
    "theta"   => 1e+6,
    "eta"     => 1e-4,
)

#=
### *Customized Structs* : *AC Solver*
=#

"""
    AbstractSolver

An abstract type representing the solver for analytic continuation
problem. It is used to build the internal type system. All the other
solvers are its subtypes.

It has the following subtypes:
* MaxEntSolver
* BarRatSolver
* NevanACSolver
* StochACSolver
* StochSKSolver
* StochOMSolver
* StochPXSolver
"""
abstract type AbstractSolver end

"""
    MaxEntSolver

It represents the analytic continuation solver that implements the
maximum entropy method.

This solver is **highly recommended**.

### Features
* Both off-diagonal and diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* The spectral weights can be negative.
* Good robustness with respect to noisy correlators.
* Numerically stable.
* Fast.

### Known Limitations
* It always tends to generate broad and smooth spectra.
* The fine features in the spectra could be smeared out.
* Sometimes it would generate keen-edged and unphysical peaks.
* When the spectra are sharp or discrete, this solver usually fails.

### Implementation

See `src/maxent.jl`.

### Tests

See `test/basic`.
"""
struct MaxEntSolver <: AbstractSolver end

"""
    BarRatSolver

It represents the analytic continuation solver that implements the
barycentric rational function approximation method.

This solver is **highly recommended**.

### Features
* Both off-diagonal and diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* Both continuous and discrete spectra are supported.
* The spectral weights can be negative.
* Moderate robustness with respect to noisy correlators.
* Numerically stable.
* It can provide analytic expressions to approximate the correlators.
* Extremely fast.

### Known Limitations
* The imaginary-time correlation correlators are not supported.
* Sometimes it would generate keen-edged and unphysical peaks.

### Implementation

See `src/rfa.jl`.

### Tests

See `test/rfa`.
"""
struct BarRatSolver <: AbstractSolver end

"""
    NevanACSolver

It represents the analytic continuation solver that implements the
Nevanlinna analytical continuation (It doesn't support the analytic
confinuations for matrix-valued Green's functions).

This solver is **not recommended**.

This solver is not well tested. If you would like to use the Nevanlinna
analytical continuation method, perhaps the official `Nevanlinna.jl`
is a better choice.

### Features
* Only diagonal Green's functions are supported.
* Only fermionic correlators are supported.
* Both continuous and discrete spectra are supported.
* It is quite accurate if the input correlators are noise-free.

### Known Limitations
* It needs Matsubara data.
* It is quite slow if Hardy basis optimization is activated.
* It doesn't support bosonic correlators directly.
* It doesn't support off-diagonal Green's functions.
* It doesn't suit the noisy quantum Monte Carlo data.
* It is numerically unstable when the input correlators are noisy.

### Implementation

See `src/nac.jl`.

### Tests

See `test/nac`.
"""
struct NevanACSolver <: AbstractSolver end

"""
    StochACSolver

It represents the analytic continuation solver that implements the
stochastic analytic continuation method (K. S. D. Beach's version).

This solver is **moderately recommended**. It is an alternative of
the `StochSK` solver.

### Features
* Only diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* Both continuous and discrete spectra are supported.
* Both imaginary-time and Matsubara data are supported.
* Good robustness with respect to noisy correlators.
* It supports the constrained sampling algorithm.
* Numerically stable.
* It is parallelized.

### Known Limitations
* It is quite slow.
* It is less accurate than the `BarRat` solver at most cases.
* It doesn't support off-diagonal Green's functions.
* It doesn't support negative spectral weights.
* If there are multiple δ-like peaks, it is not good.

### Implementation

See `src/sac.jl`.

### Tests

See `test/basic`.
"""
struct StochACSolver <: AbstractSolver end

"""
    StochSKSolver

It represents the analytic continuation solver that implements the
stochastic analytic continuation method (A. W. Sandvik's version).

This solver is **moderately recommended**. It is an alternative of
the `StochAC` solver.

### Features
* Only diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* Both continuous and discrete spectra are supported.
* Both imaginary-time and Matsubara data are supported.
* Good robustness with respect to noisy correlators.
* It supports the constrained sampling algorithm.
* Numerically stable.
* It is parallelized.

### Known Limitations
* It is quite slow.
* It is less accurate than the `BarRat` solver at most cases.
* It doesn't support off-diagonal Green's functions.
* It doesn't support negative spectral weights.
* If there are multiple δ-like peaks, it is not good.

### Implementation

See `src/san.jl`.

### Tests

See `test/basic`.
"""
struct StochSKSolver <: AbstractSolver end

"""
    StochOMSolver

It represents the analytic continuation solver that implements the
stochastic optimization method.

This solver is **less recommended**.

### Features
* Only diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* Only continuous spectra are supported.
* Good robustness with respect to noisy correlators.
* It supports the constrained sampling algorithm.
* Numerically stable.
* It is parallelized.

### Known Limitations
* It is extremely slow.
* It is less accurate than the `BarRat` solver at most cases.
* It doesn't support off-diagonal Green's functions.
* It doesn't support negative spectral weights.
* If there are multiple δ-like peaks, it is not good.

### Implementation

See `src/som.jl`.

### Tests

See `test/basic`.
"""
struct StochOMSolver <: AbstractSolver end

"""
    StochPXSolver

It represents the analytic continuation solver that implements the
stochastic pole expansion method.

This solver is **highly recommended**.

### Features
* Both off-diagonal and diagonal Green's functions are supported.
* Both fermionic and bosonic correlators are supported.
* Both continuous and discrete spectra are supported.
* The spectral weights can be negative.
* Good robustness with respect to noisy correlators.
* It supports the constrained sampling algorithm.
* Numerically stable.
* It can provide analytic expressions to approximate the correlators.
* It is parallelized.

### Known Limitations
* It is quite slow.
* It is less accurate than the `BarRat` solver at most cases.
* It requires Matsubara data as input.
* If the spectra are discrete, `a priori` knowledge is essential.
* If the spectra exhibit negative weights, `a priori` knowledge is essential.

### Implementation

See `src/spx.jl`.

### Tests

See `test/pole` and `test/lqcd`.
"""
struct StochPXSolver <: AbstractSolver end

#=
### *Customized Structs* : *Input Data*
=#

"""
    AbstractData

An abstract type representing the input data in imaginary axis. It is
used to build the internal type system.
"""
abstract type AbstractData end

"""
    RawData

Mutable struct. It represent the raw input data. The datatype `T` of raw
data may be `Float64` or `ComplexF64`.

### Members
* _grid -> Raw grid for the input data, such as τ or iωₙ.
* value -> Raw input data, such as G(τ), G(iωₙ), or Σ(iωₙ).
* error -> Error bar (standard deviation) of raw input data, σ.

See also: [`GreenData`](@ref).
"""
mutable struct RawData{T} <: AbstractData
    _grid :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

"""
    GreenData

Mutable struct. It represents the preprocessed input data. Note that it
should support arbitrary precision via `T`

### Members
* value -> Preprocessed input data.
* error -> Preprocessed error bar.
* covar -> Diagonal parts of the covariance matrix, σ².

See also: [`RawData`](@ref).
"""
mutable struct GreenData{T} <: AbstractData
    value :: Vector{T}
    error :: Vector{T}
    covar :: Vector{T}
end

#=
### *Customized Structs* : *Input Grid*
=#

"""
    AbstractGrid

An abstract type representing the imaginary axis. It is used to build
the internal type system.
"""
abstract type AbstractGrid end

"""
    FermionicImaginaryTimeGrid

Mutable struct. It represents the fermionic imaginary time grid.

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points， τᵢ.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
mutable struct FermionicImaginaryTimeGrid{T} <: AbstractGrid
    ntime :: I64
    β :: T
    τ :: Vector{T}
end

"""
    FermionicFragmentTimeGrid

Mutable struct. It represents part of the fermionic imaginary time grid.
In other words, the grid might be fragmentary。

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points, τᵢ.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
mutable struct FermionicFragmentTimeGrid{T} <: AbstractGrid
    ntime :: I64
    β :: T
    τ :: Vector{T}
end

"""
    FermionicMatsubaraGrid

Mutable struct. It represents the fermionic Matsubara frequency grid.

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
mutable struct FermionicMatsubaraGrid{T} <: AbstractGrid
    nfreq :: I64
    β :: T
    ω :: Vector{T}
end

"""
    FermionicFragmentMatsubaraGrid

Mutable struct. It represents part of the fermionic Matsubara frequency
grid. In other words, the grid might be fragmentary。

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
mutable struct FermionicFragmentMatsubaraGrid{T} <: AbstractGrid
    nfreq :: I64
    β :: T
    ω :: Vector{T}
end

"""
    BosonicImaginaryTimeGrid

Mutable struct. It represents the bosonic imaginary time grid.

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points, τᵢ.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
mutable struct BosonicImaginaryTimeGrid{T} <: AbstractGrid
    ntime :: I64
    β :: T
    τ :: Vector{T}
end

"""
    BosonicFragmentTimeGrid

Mutable struct. It represents part of the bosonic imaginary time grid.
In other words, the grid might be fragmentary。

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points, τᵢ.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
mutable struct BosonicFragmentTimeGrid{T} <: AbstractGrid
    ntime :: I64
    β :: T
    τ :: Vector{T}
end

"""
    BosonicMatsubaraGrid

Mutable struct. It represents the bosonic Matsubara frequency grid.

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
mutable struct BosonicMatsubaraGrid{T} <: AbstractGrid
    nfreq :: I64
    β :: T
    ω :: Vector{T}
end

"""
    BosonicFragmentMatsubaraGrid

Mutable struct. It represents part of the bosonic Matsubara frequency
grid. In other words, the grid might be fragmentary。 However, the first
frequency point should be present (ωₙ ≡ 0.0).

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
mutable struct BosonicFragmentMatsubaraGrid{T} <: AbstractGrid
    nfreq :: I64
    β :: T
    ω :: Vector{T}
end

#=
### *Customized Structs* : *Output Mesh*
=#

"""
    AbstractMesh

An abstract type representing the real axis. It is used to build the
internal type system.
"""
abstract type AbstractMesh end

"""
    LinearMesh

Mutable struct. A linear and uniform mesh.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`TangentMesh`](@ref).
"""
mutable struct LinearMesh{T} <: AbstractMesh
    nmesh :: I64
    wmax :: T
    wmin :: T
    mesh :: Vector{T}
    weight :: Vector{T}
end

"""
    TangentMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LinearMesh`](@ref).
"""
mutable struct TangentMesh{T} <: AbstractMesh
    nmesh :: I64
    wmax :: T
    wmin :: T
    mesh :: Vector{T}
    weight :: Vector{T}
end

"""
    LorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`HalfLorentzMesh`](@ref).
"""
mutable struct LorentzMesh{T} <: AbstractMesh
    nmesh :: I64
    wmax :: T
    wmin :: T
    mesh :: Vector{T}
    weight :: Vector{T}
end

"""
    HalfLorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on positive half-axis only.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value). It must be 0.0.
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LorentzMesh`](@ref).
"""
mutable struct HalfLorentzMesh{T} <: AbstractMesh
    nmesh :: I64
    wmax :: T
    wmin :: T
    mesh :: Vector{T}
    weight :: Vector{T}
end

"""
    DynamicMesh

Mutable struct. A mesh used internally in stochastic methods. It supports
both uniform and non-uniform meshes. The mesh is usually generated by
`util/gmesh.jl`, saved in `fmesh.inp`, and loaded dynamically during the
initialization step. This mesh should not be used as a regular mesh for
describing the spectral functions.

Note that only the StochAC and StochPX solvers support DynamicMesh. This
is because δ-like peaks in the StochAC solver and poles in the StochPX
solvers can distribute in a non-uniform mesh. While in the StochSK solver,
the δ-like peaks must be distributed in an uniform mesh, so it doesn't
support DynamicMesh. See sac.jl/calc_fmesh() and spx.jl/calc_fmesh() for
more details.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LinearMesh`](@ref).
"""
mutable struct DynamicMesh{T} <: AbstractMesh
    nmesh :: I64
    wmax :: T
    wmin :: T
    mesh :: Vector{T}
    weight :: Vector{T}
end

#=
### *Customized Structs* : *Monte Carlo Engine*
=#

"""
    AbstractMC

An abstract type representing the Monte Carlo engine. It is used to build
the internal type system.
"""
abstract type AbstractMC end

"""
    StochACMC

Mutable struct. It is used within the StochAC solver. It includes random
number generator and some counters.

### Members
* rng  -> Random number generator.
* Macc -> Counter for move operation (accepted).
* Mtry -> Counter for move operation (tried).
* Sacc -> Counter for swap operation (accepted).
* Stry -> Counter for swap operation (tried).

See also: [`StochACSolver`](@ref).
"""
mutable struct StochACMC <: AbstractMC
    rng :: AbstractRNG
    Macc :: Vector{I64}
    Mtry :: Vector{I64}
    Sacc :: Vector{I64}
    Stry :: Vector{I64}
end

"""
    StochSKMC

Mutable struct. It is used within the StochSK solver. It includes random
number generator and some counters.

### Members
* rng  -> Random number generator.
* Sacc -> Counter for single-updated operation (accepted).
* Stry -> Counter for single-updated operation (tried).
* Pacc -> Counter for pair-updated operation (accepted).
* Ptry -> Counter for pair-updated operation (tried).
* Qacc -> Counter for quadruple-updated operation (accepted).
* Qtry -> Counter for quadruple-updated operation (tried).

See also: [`StochSKSolver`](@ref).
"""
mutable struct StochSKMC
    rng :: AbstractRNG
    Sacc :: I64
    Stry :: I64
    Pacc :: I64
    Ptry :: I64
    Qacc :: I64
    Qtry :: I64
end

"""
    StochOMMC

Mutable struct. It is used within the StochOM solver. It includes random
number generator and some counters.

Because the StochOM solver supports many Monte Carlo updates, so `Macc`
and `Mtry` are vectors.

### Members
* rng  -> Random number generator.
* Macc -> Counter for move operation (accepted).
* Mtry -> Counter for move operation (tried).

See also: [`StochOMSolver`](@ref).
"""
mutable struct StochOMMC <: AbstractMC
    rng :: AbstractRNG
    Macc :: Vector{I64}
    Mtry :: Vector{I64}
end

"""
    StochPXMC

Mutable struct. It is used within the StochPX solver. It includes random
number generator and some counters.

### Members
* rng  -> Random number generator.
* Sacc -> Counter for position-updated (type 1) operation (accepted).
* Stry -> Counter for position-updated (type 1) operation (tried).
* Pacc -> Counter for position-updated (type 2) operation (accepted).
* Ptry -> Counter for position-updated (type 2) operation (tried).
* Aacc -> Counter for amplitude-updated operation (accepted).
* Atry -> Counter for amplitude-updated operation (tried).
* Xacc -> Counter for exchange operation (accepted).
* Xtry -> Counter for exchange operation (tried).

See also: [`StochPXSolver`](@ref).
"""
mutable struct StochPXMC <: AbstractMC
    rng :: AbstractRNG
    Sacc :: I64
    Stry :: I64
    Pacc :: I64
    Ptry :: I64
    Aacc :: I64
    Atry :: I64
    Xacc :: I64
    Xtry :: I64
end
