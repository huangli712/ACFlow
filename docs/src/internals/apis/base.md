# Core

*Provide basic user's interfaces for the ACFlow package.*

## Functions

```@docs
solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})
solve(grid::Vector{F64}, Gval::Vector{T}, err::T)
solve(grid::Vector{F64}, Gval::Vector{T})
solve(rd::RawData)
```

```@docs
reprod
kramers
setup_param
read_param
read_data
make_data
make_grid
make_mesh
make_model
make_kernel
```
