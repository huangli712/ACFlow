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
constraints
try_move_s
try_move_p
try_insert
try_remove
try_shift
try_width
try_height
try_split
try_merge
try_move1
try_move2
Pdx
```
