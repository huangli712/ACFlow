# Grid

*Define various grids for input data.*

## Types

```@docs
AbstractGrid
FermionicImaginaryTimeGrid
FermionicMatsubaraGrid
BosonicImaginaryTimeGrid
BosonicMatsubaraGrid
```

## Functions

### Constructors

```@docs
FermionicImaginaryTimeGrid(ntime::I64, β::F64)
FermionicMatsubaraGrid(nfreq::I64, β::F64)
```

### Base.* functions

```@docs
Base.length(fg::FermionicImaginaryTimeGrid)
Base.length(fg::FermionicMatsubaraGrid)

Base.iterate(fg::FermionicImaginaryTimeGrid)
Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)

Base.iterate(fg::FermionicMatsubaraGrid, i::I64)
Base.iterate(fg::FermionicMatsubaraGrid)


Base.eachindex(fg::FermionicImaginaryTimeGrid)
Base.eachindex(fg::FermionicMatsubaraGrid)

Base.firstindex(fg::FermionicImaginaryTimeGrid)
Base.firstindex(fg::FermionicMatsubaraGrid)

Base.lastindex(fg::FermionicImaginaryTimeGrid)
Base.lastindex(fg::FermionicMatsubaraGrid)

Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)

Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})
Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})
```

### Utilities

```@docs
rebuild
```
