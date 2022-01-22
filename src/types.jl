abstract type AbstractSolver end
struct MaxEntSolver <: AbstractSolver end
struct StochOMSolver <: AbstractSolver end
struct StochACSolver <: AbstractSolver end

abstract type AbstractModel end
struct FlatModel <: AbstractModel end
struct GaussianModel <: AbstractModel end

abstract type AbstractMesh end
struct UniformMesh <: AbstractMesh end
struct NonUniformMesh <: AbstractMesh end

abstract type AbstractGrid end
struct FermionicImaginaryTimeGrid <: AbstractGrid end
struct FermionicMatsubaraGrid <: AbstractGrid end
struct BosonicImaginaryTimeGrid <: AbstractGrid end
struct BosonicMatsubaraGrid <: AbstractGrid end

abstract type AbstractKernel end
struct FermionicImaginaryTimeKernel <: AbstractKernel end
struct FermionicMatsubaraKernel <: AbstractKernel end
struct BosonicImaginaryTimeKernel <: AbstractKernel end
struct BosonicMatsubaraKernel <: AbstractKernel end

struct RawData{T}
    mesh  :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end