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

const PMaxEnt  = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized α parameter"],
    "nalph"   => [missing, 1, :I64   , "Number of the α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we blur the kernel and spectrum"],
)

const PStochOM = Dict{String,ADT}(
)

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

abstract type AbstractData end

mutable struct RawData{T} <: AbstractData
    _grid :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

mutable struct GreenData <: AbstractData
    value :: Vector{F64}
    error :: Vector{F64}
    covar :: Vector{F64}
end

abstract type AbstractGrid end

mutable struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

mutable struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

mutable struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

mutable struct BosonicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

abstract type AbstractMesh end

mutable struct LinearMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

mutable struct TangentMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end