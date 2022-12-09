#
# Project : Gardenia
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/10
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
    "solver"  => [missing, 1, :String, "Solver for the analytical continuation problem"],
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
    "nalph"   => [missing, 1, :I64   , "Total number of the chosen α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we preblur the kernel and spectrum"],
)

# Default parameters for PMaxEnt
const _PMaxEnt = Dict{String,Any}(
    "method"  => "chi2kink",
    "nalph"   => 12,
    "alpha"   => 1e9,
    "ratio"   => 10.0,
    "blur"    => -1.0, # Negative value means off.
)

"""
    PStochAC

Dictionary for configuration parameters:
the stochastic analytical continuation method (K. S. D. Beach's version).
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
the stochastic analytical continuation method (A. W. Sandvik's version).
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
const _PStochOM = Dict{String,Any}(
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
the stochastic pole expansion.
"""
const PStochPX = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to evaluate the final spectral density"],
    "nfine"   => [missing, 1, :I64   , "Number of points of a very fine linear mesh"],
    "npole"   => [missing, 1, :I64   , "Number of poles"],
    "ntry"    => [missing, 1, :I64   , "Number of attempts (tries) to seek the solution"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo steps per attempt / try"],
    "theta"   => [missing, 1, :F64   , "Artificial inverse temperature"],
    "eta"     => [missing, 1, :F64   , "Tiny distance from the real axis"],
)

# Default parameters for PStochPX
const _PStochPX = Dict{String,Any}(
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

An abstract type representing the solver for analytical continuation
problem. It is used to build the internal type system. All the other
solvers are its sub-types.
"""
abstract type AbstractSolver end

"""
    MaxEntSolver

It represents the analytical continuation solver that implements the
maximum entropy method.
"""
struct MaxEntSolver <: AbstractSolver end

"""
    StochACSolver

It represents the analytical continuation solver that implements the
stochastic analytical continuation method (K. S. D. Beach's version).
"""
struct StochACSolver <: AbstractSolver end

"""
    StochSKSolver

It represents the analytical continuation solver that implements the
stochastic analytical continuation method (A. W. Sandvik's version).
"""
struct StochSKSolver <: AbstractSolver end

"""
    StochOMSolver

It represents the analytical continuation solver that implements the
stochastic optimization method.
"""
struct StochOMSolver <: AbstractSolver end

"""
    StochPXSolver

It represents the analytical continuation solver that implements the
stochastic pole expansion.
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

Mutable struct. It represent the raw input data. The datatype of raw data
may be float or complex.

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

Mutable struct. It represents the preprocessed input data.

### Members

* value -> Preprocessed input data.
* error -> Preprocessed error bar.
* covar -> Diagonal parts of the covariance matrix, σ².

See also: [`RawData`](@ref).
"""
mutable struct GreenData <: AbstractData
    value :: Vector{F64}
    error :: Vector{F64}
    covar :: Vector{F64}
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
* τ     -> Vector of grid points.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
mutable struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    FermionicMatsubaraGrid

Mutable struct. It represents the fermionic Matsubara frequency grid.

### Members

* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
mutable struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

"""
    BosonicImaginaryTimeGrid

Mutable struct. It represents the bosonic imaginary time grid.

### Members

* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
mutable struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    BosonicMatsubaraGrid

Mutable struct. It represents the bosonic Matsubara frequency grid.

### Members

* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
mutable struct BosonicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
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

* nmesh  -> Number of mesh points
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`TangentMesh`](@ref).
"""
mutable struct LinearMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    TangentMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members

* nmesh  -> Number of mesh points
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LinearMesh`](@ref).
"""
mutable struct TangentMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    LorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members

* nmesh  -> Number of mesh points
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`HalfLorentzMesh`](@ref).
"""
mutable struct LorentzMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    HalfLorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on positive half-axis only.

### Members

* nmesh  -> Number of mesh points
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value). It must be 0.0.
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LorentzMesh`](@ref).
"""
mutable struct HalfLorentzMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
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
