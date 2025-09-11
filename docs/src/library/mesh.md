*Meshes on real axis.*

The spectral functions are always defined on real axis. The ACFlow toolkit supports various uniform and non-uniform meshes. In order to build these meshes, we need some additional control parameters, including ``f_1`` and `cut`. They should be setup by using the parameter [`pmesh`](@ref pmesh).

```@index
Pages = ["mesh.md"]
```

## Types

```@docs
AbstractMesh
LinearMesh
TangentMesh
LorentzMesh
HalfLorentzMesh
DynamicMesh
```

## Constructors

```@docs
LinearMesh(nmesh::I64, wmin::F64, wmax::F64)
TangentMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 2.1)
LorentzMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 0.01)
HalfLorentzMesh(nmesh::I64, wmax::F64, 𝑝::F64 = 0.01)
DynamicMesh(mesh::Vector{F64})
```

## Base.* Functions

```@docs
Base.length(am::AbstractMesh)
Base.iterate(am::AbstractMesh)
Base.iterate(am::AbstractMesh, i::I64)
Base.eachindex(am::AbstractMesh)
Base.firstindex(am::AbstractMesh)
Base.lastindex(am::AbstractMesh)
Base.getindex(am::AbstractMesh, ind::I64)
Base.getindex(am::AbstractMesh, I::UnitRange{I64})
```

## Utilities

```@docs
nearest
```
