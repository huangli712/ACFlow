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

```@docs
LinearMesh(nmesh::I64, wmin::F64, wmax::F64)
Base.length
Base.iterate
Base.getindex
Base.eachindex
Base.firstindex
Base.lastindex
nearest
```
