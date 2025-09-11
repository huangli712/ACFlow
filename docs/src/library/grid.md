*Grids on imaginary axis.*

In the ACFlow toolkit, the input correlators are defined on `grid`, while the calculated spectral functions are defined on `mesh`. The ACFlow toolkit supports both imaginary time and Matsubara frequency grids. Note that for Matsubara frequency grid, the bosonic and fermionic grids are different.

```@index
Pages = ["grid.md"]
```

## Types

```@docs
AbstractGrid
FermionicImaginaryTimeGrid
FermionicFragmentTimeGrid
FermionicMatsubaraGrid
FermionicFragmentMatsubaraGrid
BosonicImaginaryTimeGrid
BosonicFragmentTimeGrid
BosonicMatsubaraGrid
BosonicFragmentMatsubaraGrid
```

## Constructors

```@docs
FermionicImaginaryTimeGrid(ntime::I64, β::F64)
FermionicFragmentTimeGrid(β::F64, τ::Vector{F64})
FermionicMatsubaraGrid(nfreq::I64, β::F64)
FermionicFragmentMatsubaraGrid(β::F64, ω::Vector{F64})
BosonicImaginaryTimeGrid(ntime::I64, β::F64)
BosonicFragmentTimeGrid(β::F64, τ::Vector{F64})
BosonicMatsubaraGrid(nfreq::I64, β::F64)
BosonicFragmentMatsubaraGrid(β::F64, ω::Vector{F64})
```

## Base.* Functions

```@docs
Base.length(fg::FermionicImaginaryTimeGrid)
Base.length(fg::FermionicFragmentTimeGrid)
Base.length(fg::FermionicMatsubaraGrid)
Base.length(fg::FermionicFragmentMatsubaraGrid)
Base.length(bg::BosonicImaginaryTimeGrid)
Base.length(bg::BosonicFragmentTimeGrid)
Base.length(bg::BosonicMatsubaraGrid)
Base.length(bg::BosonicFragmentMatsubaraGrid)
Base.iterate(fg::FermionicImaginaryTimeGrid)
Base.iterate(fg::FermionicFragmentTimeGrid)
Base.iterate(fg::FermionicMatsubaraGrid)
Base.iterate(fg::FermionicFragmentMatsubaraGrid)
Base.iterate(bg::BosonicImaginaryTimeGrid)
Base.iterate(bg::BosonicFragmentTimeGrid)
Base.iterate(bg::BosonicMatsubaraGrid)
Base.iterate(bg::BosonicFragmentMatsubaraGrid)
Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)
Base.iterate(fg::FermionicFragmentTimeGrid, i::I64)
Base.iterate(fg::FermionicMatsubaraGrid, i::I64)
Base.iterate(fg::FermionicFragmentMatsubaraGrid, i::I64)
Base.iterate(bg::BosonicImaginaryTimeGrid, i::I64)
Base.iterate(bg::BosonicFragmentTimeGrid, i::I64)
Base.iterate(bg::BosonicMatsubaraGrid, i::I64)
Base.iterate(bg::BosonicFragmentMatsubaraGrid, i::I64)
Base.eachindex(fg::FermionicImaginaryTimeGrid)
Base.eachindex(fg::FermionicFragmentTimeGrid)
Base.eachindex(fg::FermionicMatsubaraGrid)
Base.eachindex(fg::FermionicFragmentMatsubaraGrid)
Base.eachindex(bg::BosonicImaginaryTimeGrid)
Base.eachindex(bg::BosonicFragmentTimeGrid)
Base.eachindex(bg::BosonicMatsubaraGrid)
Base.eachindex(bg::BosonicFragmentMatsubaraGrid)
Base.firstindex(fg::FermionicImaginaryTimeGrid)
Base.firstindex(fg::FermionicFragmentTimeGrid)
Base.firstindex(fg::FermionicMatsubaraGrid)
Base.firstindex(fg::FermionicFragmentMatsubaraGrid)
Base.firstindex(bg::BosonicImaginaryTimeGrid)
Base.firstindex(bg::BosonicFragmentTimeGrid)
Base.firstindex(bg::BosonicMatsubaraGrid)
Base.firstindex(bg::BosonicFragmentMatsubaraGrid)
Base.lastindex(fg::FermionicImaginaryTimeGrid)
Base.lastindex(fg::FermionicFragmentTimeGrid)
Base.lastindex(fg::FermionicMatsubaraGrid)
Base.lastindex(fg::FermionicFragmentMatsubaraGrid)
Base.lastindex(bg::BosonicImaginaryTimeGrid)
Base.lastindex(bg::BosonicFragmentTimeGrid)
Base.lastindex(bg::BosonicMatsubaraGrid)
Base.lastindex(bg::BosonicFragmentMatsubaraGrid)
Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
Base.getindex(fg::FermionicFragmentTimeGrid, ind::I64)
Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
Base.getindex(fg::FermionicFragmentMatsubaraGrid, ind::I64)
Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
Base.getindex(bg::BosonicFragmentTimeGrid, ind::I64)
Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
Base.getindex(bg::BosonicFragmentMatsubaraGrid, ind::I64)
Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})
Base.getindex(fg::FermionicFragmentTimeGrid, I::UnitRange{I64})
Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})
Base.getindex(fg::FermionicFragmentMatsubaraGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicImaginaryTimeGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicFragmentTimeGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicMatsubaraGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicFragmentMatsubaraGrid, I::UnitRange{I64})
```

## Utilities

```@docs
rebuild!
resize!
reverse!
```
