# Grids

*Grids on imaginary axis.*

In the **ACFlow** package, the input correlators are defined on `grid`, while the calculated spectral functions are defined on `mesh`. The **ACFlow** package supports both imaginary-time and Matsubara frequency grids.

## Contents

```@contents
Pages = ["grid.md"]
Depth = 2
```

## Index

```@index
Pages = ["grid.md"]
```

## Types

```@docs
AbstractGrid
FermionicImaginaryTimeGrid
FermionicMatsubaraGrid
BosonicImaginaryTimeGrid
BosonicMatsubaraGrid
```

## Constructors

```@docs
FermionicImaginaryTimeGrid(ntime::I64, β::F64)
FermionicMatsubaraGrid(nfreq::I64, β::F64)
BosonicImaginaryTimeGrid(ntime::I64, β::F64)
BosonicMatsubaraGrid(nfreq::I64, β::F64)
```

## Base.* Functions

```@docs
Base.length(fg::FermionicImaginaryTimeGrid)
Base.length(fg::FermionicMatsubaraGrid)
Base.length(bg::BosonicImaginaryTimeGrid)
Base.length(bg::BosonicMatsubaraGrid)
Base.iterate(fg::FermionicImaginaryTimeGrid)
Base.iterate(fg::FermionicMatsubaraGrid)
Base.iterate(bg::BosonicImaginaryTimeGrid)
Base.iterate(bg::BosonicMatsubaraGrid)
Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)
Base.iterate(fg::FermionicMatsubaraGrid, i::I64)
Base.iterate(bg::BosonicImaginaryTimeGrid, i::I64)
Base.iterate(bg::BosonicMatsubaraGrid, i::I64)
Base.eachindex(fg::FermionicImaginaryTimeGrid)
Base.eachindex(fg::FermionicMatsubaraGrid)
Base.eachindex(bg::BosonicImaginaryTimeGrid)
Base.eachindex(bg::BosonicMatsubaraGrid)
Base.firstindex(fg::FermionicImaginaryTimeGrid)
Base.firstindex(fg::FermionicMatsubaraGrid)
Base.firstindex(bg::BosonicImaginaryTimeGrid)
Base.firstindex(bg::BosonicMatsubaraGrid)
Base.lastindex(fg::FermionicImaginaryTimeGrid)
Base.lastindex(fg::FermionicMatsubaraGrid)
Base.lastindex(bg::BosonicImaginaryTimeGrid)
Base.lastindex(bg::BosonicMatsubaraGrid)
Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})
Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicImaginaryTimeGrid, I::UnitRange{I64})
Base.getindex(bg::BosonicMatsubaraGrid, I::UnitRange{I64})
```

## Utilities

```@docs
rebuild
```
