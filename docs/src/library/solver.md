# Solvers

*Define various solvers for the ACFlow toolkit.*

Now the ACFlow toolkit supports five analytical continuation solvers. They are:

* `MaxEnt` (Maximum Entropy Method, see `maxent.jl`)
* `StochAC` (Stochastic Analytical Continuations, see `sac.jl`)
* `StochSK` (Stochastic Analytical Continuation, see `san.jl`)
* `StochOM` (Stochastic Optimization Method, see `som.jl`)
* `StochPX` (Stochastic Pole Expansion, see `spx.jl`)

Note that the `StochAC` solver is based on the Beach's variant, while the `StochSK` solver is based on the Sandvik's variant.

## Contents

```@contents
Pages = ["solver.md"]
Depth = 3
```

## Index

```@index
Pages = ["solver.md"]
```

## Abstract Structs

```@docs
AbstractSolver
AbstractMC
```

## MaxEnt Solver

### Structs

```@docs
MaxEntSolver
MaxEntContext
```

### Functions

```@docs
solve(S::MaxEntSolver, rd::RawData)
init(S::MaxEntSolver, rd::RawData)
ACFlow.run(mec::MaxEntContext)
ACFlow.last(mec::MaxEntContext, svec::Vector, sol::Dict)
historic
classic
bryan
chi2kink
optimizer
precompute
f_and_J
f_and_J_offdiag
svd_to_real
svd_to_real_offdiag
calc_entropy
calc_entropy_offdiag
calc_bayes
calc_bayes_offdiag
calc_chi2(mec::MaxEntContext, A::Vector{F64})
```

## StochAC Solver

### Structs

```@docs
StochACSolver
StochACMC
StochACElement
StochACContext
```

### Functions

```@docs
solve(S::StochACSolver, rd::RawData)
init(S::StochACSolver, rd::RawData)
ACFlow.run(MC::StochACMC, SE::StochACElement, SC::StochACContext)
prun(S::StochACSolver, p1::Dict{String,Vector{Any}}, p2::Dict{String,Vector{Any}}, MC::StochACMC, SE::StochACElement, SC::StochACContext)
average(step::F64, SC::StochACContext)
ACFlow.last(SC::StochACContext, Aout::Array{F64,2}, Uα::Vector{F64})
warmup(MC::StochACMC, SE::StochACElement, SC::StochACContext)
sample(MC::StochACMC, SE::StochACElement, SC::StochACContext)
measure(SE::StochACElement, SC::StochACContext)
init_mc(S::StochACSolver)
init_element(S::StochACSolver, rng::AbstractRNG, allow::Vector{I64})
init_iodata(S::StochACSolver, rd::RawData)
calc_fmesh(S::StochACSolver)
calc_phi
calc_delta
calc_hamil
calc_htau
calc_alpha
constraints(S::StochACSolver)
try_move_a(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
try_move_p(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
try_move_x(MC::StochACMC, SE::StochACElement, SC::StochACContext)
```

## StochSK Solver

### Structs

```@docs
StochSKSolver
StochSKMC
StochSKElement
StochSKContext
```

### Functions

```@docs
solve(S::StochSKSolver, rd::RawData)
init(S::StochSKSolver, rd::RawData)
ACFlow.run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
prun(S::StochSKSolver, p1::Dict{String,Vector{Any}}, p2::Dict{String,Vector{Any}}, MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
average(step::F64, SC::StochSKContext)
ACFlow.last(SC::StochSKContext, Asum::Vector{F64}, χ²vec::Vector{F64}, Θvec::Vector{F64})
warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
measure(SE::StochSKElement, SC::StochSKContext)
ACFlow.shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
init_mc(S::StochSKSolver)
init_element(S::StochSKSolver, rng::AbstractRNG, allow::Vector{I64})
init_iodata(S::StochSKSolver, rd::RawData)
calc_fmesh(S::StochSKSolver)
calc_correlator
calc_goodness
calc_theta
constraints(S::StochSKSolver)
try_move_s(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
try_move_p(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
try_move_q(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
```

## StochOM Solver

### Structs

```@docs
StochOMSolver
StochOMMC
Box
StochOMElement
StochOMContext
```

### Functions

```@docs
solve(S::StochOMSolver, rd::RawData)
init(S::StochOMSolver, rd::RawData)
ACFlow.run(MC::StochOMMC, SC::StochOMContext)
prun(S::StochOMSolver, p1::Dict{String,Vector{Any}}, p2::Dict{String,Vector{Any}}, MC::StochOMMC, SC::StochOMContext)
average(SC::StochOMContext)
ACFlow.last(SC::StochOMContext, Aout::Vector{F64})
update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)
init_mc(S::StochOMSolver)
init_element(MC::StochOMMC, SC::StochOMContext)
init_iodata(S::StochOMSolver, rd::RawData)
init_context(S::StochOMSolver)
calc_lambda
calc_error
calc_green(Λ::Array{F64,2}, nk::I64)
calc_norm
constraints(e₁::F64, e₂::F64)
try_insert
try_remove
try_shift
try_width
try_height
try_split
try_merge
Pdx
```

## StochPX Solver

### Structs

```@docs
StochPXSolver
StochPXMC
StochPXElement
StochPXContext
```

### Functions

```@docs
solve(S::StochPXSolver, rd::RawData)
init(S::StochPXSolver, rd::RawData)
ACFlow.run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
prun(S::StochPXSolver, p1::Dict{String,Vector{Any}}, p2::Dict{String,Vector{Any}}, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
average(SC::StochPXContext)
ACFlow.last(SC::StochPXContext, Aout::Vector{F64}, Gout::Vector{C64}, Gᵣ::Vector{F64})
sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
measure(t::I64, SE::StochPXElement, SC::StochPXContext)
init_mc(S::StochPXSolver)
init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})
init_iodata(S::StochPXSolver, rd::RawData)
init_context(S::StochPXSolver)
reset_mc(MC::StochPXMC)
reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
calc_fmesh(S::StochPXSolver)
calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh)
calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
calc_green(P::Vector{I64}, A::Vector{F64}, Λ::Array{F64,2})
calc_green(P::Vector{I64}, A::Vector{F64}, mesh::AbstractMesh, fmesh::AbstractMesh)
calc_green(P::Vector{I64}, A::Vector{F64}, mesh::AbstractMesh, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
constraints(S::StochPXSolver)
try_move_s(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_x(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
```
