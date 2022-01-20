
abstract type AbstractKernel end
struct FermionicImaginaryTimeKernel <: AbstractKernel end
struct FermionicMatsubaraKernel <: AbstractKernel end
struct BosonicImaginaryTimeKernel <: AbstractKernel end
struct BosonicMatsubaraKernel <: AbstractKernel end