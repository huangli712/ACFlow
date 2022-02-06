#
# Project : Gardenia
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/06
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
    PCOMM
"""
const PCOMM    = Dict{String,ADT}(
    "finput"  => [missing, 1, :String, "Filename for input correlation function"],
    "solver"  => [missing, 1, :String, "Solve for the analytical continuation problem"],
    "ktype"   => [missing, 1, :String, "Type of kernel function"],
    "mtype"   => [missing, 1, :String, "Type of default model"],
    "grid"    => [missing, 1, :String, "Grid for input correlation function"],
    "mesh"    => [missing, 1, :String, "Mesh for output spectrum"],
    "ngrid"   => [missing, 1, :I64   , "Number of grid points"],
    "nmesh"   => [missing, 1, :I64   , "Number of mesh points"],
    "wmax"    => [missing, 1, :F64   , "Maximum value of mesh"],
    "wmin"    => [missing, 1, :F64   , "Minimum value of mesh"],
    "beta"    => [missing, 1, :F64   , "Inverse temperature"],
    "offdiag" => [missing, 1, :Bool  , "Is it an offdiagonal element in matrix-valued function"],
)

"""
    PMaxEnt
"""
const PMaxEnt  = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized α parameter"],
    "nalph"   => [missing, 1, :I64   , "Number of the α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we blur the kernel and spectrum"],
)

"""
    PStochAC
"""
const PStochAC = Dict{String,ADT}(
    "nfine"   => [missing, 1, :I64   , "Number of points for very fine mesh"],
    "ngamm"   => [missing, 1, :I64   , "Number of δ functions"],
    "nalph"   => [missing, 1, :I64   , "Number of the α parameters"],
    "nwarm"   => [missing, 1, :I64   , "Number of Monte Carlo warming-up steps"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo sweeping steps"],
    "ndump"   => [missing, 1, :I64   , "Intervals for monitoring Monte Carlo sweeps"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
)

"""
    PStochOM
"""
const PStochOM = Dict{String,ADT}(
)

#=
### *Customized Structs* : *AC Solver*
=#

"""
    AbstractSolver
"""
abstract type AbstractSolver end

"""
    MaxEntSolver
"""
struct MaxEntSolver <: AbstractSolver end

"""
    StochACSolver
"""
struct StochACSolver <: AbstractSolver end

"""
    StochOMSolver
"""
struct StochOMSolver <: AbstractSolver end

#=
### *Customized Structs* : *Input Data*
=#

"""
    AbstractData
"""
abstract type AbstractData end

"""
    RawData
"""
mutable struct RawData{T} <: AbstractData
    _grid :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

"""
    GreenData
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
"""
abstract type AbstractGrid end

"""
    FermionicImaginaryTimeGrid
"""
mutable struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    FermionicMatsubaraGrid
"""
mutable struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

"""
    BosonicImaginaryTimeGrid
"""
mutable struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    BosonicMatsubaraGrid
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
"""
abstract type AbstractMesh end

"""
    LinearMesh
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
"""
mutable struct TangentMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end
