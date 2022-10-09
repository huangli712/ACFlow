# Mesh

*Define various meshes for output data.*

## Types

```@docs
AbstractMesh
LinearMesh
TangentMesh
LorentzMesh
HalfLorentzMesh
```

## Functions

### Constructors

```@docs
LinearMesh(nmesh::I64, wmin::F64, wmax::F64)
LinearMesh(mesh::Vector{F64})
TangentMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 2.1)
LorentzMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 0.01)
HalfLorentzMesh(nmesh::I64, wmax::F64, 𝑝::F64 = 0.01)
```

### Base.* functions

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

### Utilities

```@docs
nearest
```
