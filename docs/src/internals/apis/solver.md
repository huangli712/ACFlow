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
run(mec::MaxEntContext)
last(mec::MaxEntContext, svec::Vector, sol::Dict)
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
