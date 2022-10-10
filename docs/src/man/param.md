# Parameters

## Contents

```@contents
Pages = ["param.md"]
Depth = 2
```

## [BASE] block

### finput

*Definition:*

Filename for the input data.

*Type*

String.

*Comment*

### solver

*Definition:*

This parameter specifies solver for the analytical continuation problem. Now **ACFlow.jl** supports four solvers. They are:

* MaxEnt
* StochAC
* StochSK
* StochOM

Here, **MaxEnt** means the maximum entropy method, **StochAC** means the stochastic analytical continuation method (K. S. D. Beach's version), **StochSK** means the stochastic analytical continuation method (A. W. Sandvik's version), and **StochOM** means the stochastic optimization method.

*Type*

String.

*Comment*

### ktype

*Definition:*

It denotes type of kernel function. Now **ACFlow.jl** supports three types of kernel functions. They are:

* fermi
* boson
* bsymm


*Type*

*Comment*

### mtype

*Definition:*

*Type*





*Comment*

### grid

*Definition:*

*Type*





*Comment*

### mesh

*Definition:*

*Type*





*Comment*

### ngrid

*Definition:*

*Type*





*Comment*

### nmesh

*Definition:*

*Type*





*Comment*

### wmax

*Definition:*

*Type*





*Comment*

### wmin

*Definition:*

*Type*





*Comment*

### beta

*Definition:*

*Type*





*Comment*

### offdiag

*Definition:*

*Type*





*Comment*

### pmodel

*Definition:*

*Type*





*Comment*

### pmesh

*Definition:*

*Type*





*Comment*

### exclude

*Definition:*

*Type*





*Comment*

## [MaxEnt]

### method

*Definition:*

*Type*





*Comment*

### nalph

*Definition:*

*Type*





*Comment*

### alpha

*Definition:*

*Type*





*Comment*

### ratio

*Definition:*

*Type*





*Comment*

### blur

*Definition:*

*Type*





*Comment*

## [StochAC] block

### nfine

*Definition:*

*Type*





*Comment*

### ngamm

*Definition:*

*Type*





*Comment*

### nwarm

*Definition:*

*Type*





*Comment*

### nstep

*Definition:*

*Type*





*Comment*

### ndump

*Definition:*

*Type*





*Comment*

### nalph

*Definition:*

*Type*





*Comment*

### alpha

*Definition:*

*Type*





*Comment*

### ratio

*Definition:*

*Type*





*Comment*

## [StochSK] block

### nfine

*Definition:*

*Type*





*Comment*

### ngamm

*Definition:*

*Type*





*Comment*

### nwarm

*Definition:*

*Type*





*Comment*

### nstep

*Definition:*

*Type*





*Comment*

### ndump

*Definition:*

*Type*





*Comment*

### retry

*Definition:*

*Type*





*Comment*

### theta

*Definition:*

*Type*





*Comment*

### ratio

*Definition:*

*Type*





*Comment*

## [StochOM] block

### ntry

*Definition:*

*Type*





*Comment*

### nstep

*Definition:*

*Type*





*Comment*

### nbox

*Definition:*

*Type*





*Comment*

### sbox

*Definition:*

*Type*





*Comment*

### wbox

*Definition:*

*Type*





*Comment*

### norm

*Definition:*

*Type*





*Comment*
