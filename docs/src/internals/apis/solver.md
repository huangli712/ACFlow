# Solvers

*Define solvers for the ACFlow package.*

## Common

### Structs

```@docs
AbstractSolver
AbstractMC
```

## MaxEntSolver

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
calc_chi2
```

## StochAC solver

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
calc_xmesh
calc_phi
calc_delta
calc_hamil
calc_htau
calc_alpha
constraints(S::StochACSolver)
try_mov1
try_mov2
try_swap
```

## StochSK solver

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
ACFlow.last(SC::StochSKContext, Asum::Vector{F64})
warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
measure(SE::StochSKElement, SC::StochSKContext)
shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
init_mc(S::StochSKSolver)
init_element(S::StochSKSolver, rng::AbstractRNG, allow::Vector{I64})
init_iodata(S::StochSKSolver, rd::RawData)
calc_fmesh(S::StochSKSolver)
calc_correlator
calc_goodness
constraints(S::StochSKSolver)
try_move_s
try_move_p
```

## StochOM solver

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
```
