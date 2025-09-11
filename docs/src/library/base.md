*Provide basic user's interfaces for the ACFlow toolkit.*

```@index
Pages = ["base.md"]
```

## Solvers

```@docs
solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T}) where {T}
solve(grid::Vector{F64}, Gval::Vector{T}, err::T) where {T}
solve(grid::Vector{F64}, Gval::Vector{T}) where {T}
solve(rd::RawData)
```

## Parameters

```@docs
setup_param
read_param
```

## Data

```@docs
read_data
make_data
```

## Grids

```@docs
make_grid
```

## Meshes

```@docs
make_mesh
```

## Models

```@docs
make_model
```

## Kernels

```@docs
make_kernel
```

## Postprocessing

```@docs
reprod
kramers
```
