abstract type AbstractMesh end
struct UniformMesh <: AbstractMesh end
struct NonUniformMesh <: AbstractMesh end

abstract type AbstractGrid end
struct FermionicImaginaryTimeGrid <: AbstractGrid end
struct FermionicFrequencyGrid <: AbstractGrid end
struct BosonicImaginaryTimeGrid <: AbstractGrid end
struct BosonicFrequencyGrid <: AbstractGrid end