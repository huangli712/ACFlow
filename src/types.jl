
abstract type AbstractData end
struct GreenData <: AbstractData
    data
    error
end

struct SigmaData <: AbstractData end
struct ChiData <: AbstractData end

abstract type AbstractGrid end
struct RealFrequencyGrid <: AbstractGrid end
struct ImaginaryTimeGrid <: AbstractGrid end
struct FermionicMatsubaraGrid <: AbstractGrid
    β
    grid
end
struct BosonicMatsubaraGrid <: AbstractGrid end
