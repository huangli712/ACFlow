
#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

const PCOMM    = Dict{String,ADT}(
    "solver" => [missing, 1, :String, ""],
    "kernel" => [missing, 1, :String, ""],
    "model"  => [missing, 1, :String, ""],
    "mesh"   => [missing, 1, :String, ""],
    "wmax"   => [missing, 1, :F64, ""],
    "wmin"   => [missing, 1, :F64, ""],
    "beta"   => [missing, 1, :F64, ""],
)

const PMaxEnt  = Dict{String,ADT}(
    "method" => [missing, 1, :String, ""],
    "alpha"  => [missing, 1, :F64, ""],
)

const PStochOM = Dict{String,ADT}(
)

const PStochAC = Dict{String,ADT}(
)

abstract type AbstractSolver end
struct MaxEntSolver <: AbstractSolver end
struct StochOMSolver <: AbstractSolver end
struct StochACSolver <: AbstractSolver end

abstract type AbstractKernel end
struct FermionicImaginaryTimeKernel <: AbstractKernel end
struct FermionicFrequencyKernel <: AbstractKernel end
struct BosonicImaginaryTimeKernel <: AbstractKernel end
struct BosonicFrequencyKernel <: AbstractKernel end

abstract type AbstractModel end
struct FlatModel <: AbstractModel end
struct GaussianModel <: AbstractModel end

abstract type AbstractMesh end
struct UniformMesh <: AbstractMesh end
struct NonUniformMesh <: AbstractMesh end

abstract type AbstractGrid end
struct FermionicImaginaryTimeGrid <: AbstractGrid end
struct FermionicFrequencyGrid <: AbstractGrid end
struct BosonicImaginaryTimeGrid <: AbstractGrid end
struct BosonicFrequencyGrid <: AbstractGrid end
