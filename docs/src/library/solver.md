*Define various solvers for the ACFlow toolkit.*

Now the ACFlow toolkit supports seven analytic continuation solvers. They are:

* `MaxEnt` (Maximum Entropy Method, see `maxent.jl`)
* `BarRat` (Barycentric Rational Function Approximation, see `rfa.jl`)
* `NevanAC` (Nevanlinna Analytical Continuation, see `nac.jl`)
* `StochAC` (Stochastic Analytic Continuation, see `sac.jl`)
* `StochSK` (Stochastic Analytic Continuation, see `san.jl`)
* `StochOM` (Stochastic Optimization Method, see `som.jl`)
* `StochPX` (Stochastic Pole Expansion, see `spx.jl`)

!!! note

    The `StochAC` solver is based on the Beach's variant, while the `StochSK` solver is based on the Sandvik's variant.

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
precompute(Gᵥ::Vector{F64}, σ²::Vector{F64}, am::AbstractMesh, D::Vector{F64}, K::Matrix{F64})
f_and_J
f_and_J_od
svd_to_real
svd_to_real_od
calc_entropy
calc_entropy_od
calc_bayes
calc_bayes_od
calc_chi2(mec::MaxEntContext, A::Vector{F64})
```

## BarRat Solver

### Structs

```@docs
BarRatSolver
BarRatContext
BarycentricFunction
PronyApproximation
```

### Functions

```@docs
solve(S::BarRatSolver, rd::RawData)
init(S::BarRatSolver, rd::RawData)
ACFlow.run(brc::BarRatContext)
ACFlow.last(brc::BarRatContext)
aaa
poles!
bc_nodes
bc_values
bc_weights
bc_degree
bc_poles
prony_data
prony_svd
prony_idx
prony_v
prony_gamma
prony_omega
```

## NevanAC Solver

### Structs

```@docs
NevanACSolver
NevanACContext
```

### Functions

```@docs
solve(S::NevanACSolver, rd::RawData)
init(S::NevanACSolver, rd::RawData)
ACFlow.run(nac::NevanACContext)
ACFlow.last(nac::NevanACContext)
precompute(grid::AbstractGrid, mesh::AbstractMesh, Gᵥ::Vector{APC})
calc_mobius
calc_inv_mobius
calc_pick
calc_phis
calc_abcd
calc_hbasis
calc_hmatrix
calc_theta(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
calc_green(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
calc_noptim
calc_hmin!
calc_hopt!
hardy_optimize!
smooth_norm
check_pick
check_causality
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
init_iodata(S::StochACSolver, rd::RawData)
init_mc(S::StochACSolver)
init_element(S::StochACSolver, rng::AbstractRNG, allow::Vector{I64})
init_context(SE::StochACElement, Gᵥ::Vector{F64}, σ¹::Vector{F64}, allow::Vector{I64}, grid::AbstractGrid, mesh::AbstractMesh, fmesh::AbstractMesh)
calc_fmesh(S::StochACSolver)
calc_phi
calc_delta
calc_hamil
calc_htau
calc_alpha
constraints(S::StochACSolver, fmesh::AbstractMesh)
try_move_s(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
try_move_p(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
try_move_a(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
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
init_iodata(S::StochSKSolver, rd::RawData)
init_mc(S::StochSKSolver)
init_element(S::StochSKSolver, rng::AbstractRNG, allow::Vector{I64})
init_context(SE::StochSKElement, Gᵥ::Vector{F64}, σ¹::Vector{F64}, allow::Vector{I64}, grid::AbstractGrid, mesh::AbstractMesh, fmesh::AbstractMesh)
calc_fmesh(S::StochSKSolver)
calc_correlator
calc_goodness
calc_theta
constraints(S::StochSKSolver, fmesh::AbstractMesh)
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
init_iodata(S::StochOMSolver, rd::RawData)
init_mc(S::StochOMSolver)
init_element(MC::StochOMMC, SC::StochOMContext)
init_context(S::StochOMSolver, grid::AbstractGrid)
eval_lambda
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
init_iodata(S::StochPXSolver, rd::RawData)
init_mc(S::StochPXSolver)
init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})
init_context(S::StochPXSolver, SE::StochPXElement, grid::AbstractGrid, fmesh::AbstractMesh, Gᵥ::Vector{F64})
reset_mc(MC::StochPXMC)
reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
calc_fmesh(S::StochPXSolver)
calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh, Gᵥ::Vector{F64})
calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh)
calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
calc_green(t::I64, SC::StochPXContext, real_axis::Bool)
calc_green(P::Vector{I64}, A::Vector{F64}, 𝕊::Vector{F64}, Λ::Array{F64,2})
calc_green(P::Vector{I64}, A::Vector{F64}, 𝕊::Vector{F64}, mesh::AbstractMesh, fmesh::AbstractMesh)
calc_green(P::Vector{I64}, A::Vector{F64}, 𝕊::Vector{F64}, mesh::AbstractMesh, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
constraints(S::StochPXSolver, fmesh::AbstractMesh)
try_move_s(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
try_move_x(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
```
