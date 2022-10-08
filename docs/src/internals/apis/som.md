# Solver

*Define solvers for the ACFlow package.*

*Source: som.jl*

## Contents

```@contents
Pages = ["som.md"]
```

## Index

```@index
Pages = ["som.md"]
```

## Structs

```@docs
Box
MaxEntContext
StochACElement
StochACContext
StochSKElement
StochSKContext
StochOMElement
StochOMContext
```

## Functions

```@docs
init
ACFlow.run
prun
average
ACFlow.last
update
warmup
sample
measure
shuffle
init_mc
init_element
init_context
init_iodata
calc_fmesh
calc_xmesh
calc_correlator
calc_goodness
calc_lambda
calc_error
calc_green
calc_norm
calc_alpha
calc_htau
calc_chi2
calc_hamil
calc_phi
calc_delta
calc_entropy
calc_entropy_offdiag
calc_bayes
calc_bayes_offdiag
classic
historic
bryan
chi2kink
constraints
f_and_J
f_and_J_offdiag
svd_to_real
svd_to_real_offdiag
precompute
try_move_s
try_move_p
try_insert
try_remove
try_shift
try_width
try_height
try_split
try_merge
try_mov1
try_mov2
try_swap
Pdx
```
