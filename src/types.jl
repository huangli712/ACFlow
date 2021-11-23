
abstract type AbstractData end
struct GreenData <: AbstractData end
struct SigmaData <: AbstractData end
struct ChiData <: AbstractData end
struct SpectrumData <: AbstractData end

abstract type AbstractStatistics end
struct FermiStatistics <: AbstractStatistics end
struct BosonStatistics <: AbstractStatistics end

abstract type AbstractGrid end
struct RealFrequencyGrid <: AbstractGrid end
struct ImaginaryTimeGrid <: AbstractGrid end
struct FermionicMatsubaraGrid <: AbstractGrid end
struct BosonicMatsubaraGrid <: AbstractGrid end
