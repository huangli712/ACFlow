# Core

*Provide basic user's interfaces for the ACFlow package.*

## Functions

### Solve analytical continuation problems

```@docs
solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})
solve(grid::Vector{F64}, Gval::Vector{T}, err::T)
solve(grid::Vector{F64}, Gval::Vector{T})
solve(rd::RawData)
```

### Parameters

```@docs
setup_param
read_param
```

### Data

```@docs
read_data
make_data
```

### Grid and mesh

```@docs
make_grid
make_mesh
```

### Models

```@docs
make_model
```

### Kernels

```@docs
make_kernel
```

### Postprocess

```@docs
reprod
kramers
```
