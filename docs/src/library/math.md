# Math

*Define some essential mathematical functions.*

## Contents

```@contents
Pages = ["math.md"]
Depth = 3
```

## Index

```@index
Pages = ["math.md"]
```

## Root finding

```@docs
secant
newton
```

## Numerical integrations

```@docs
trapz
simpson
```

## Interpolations

### Structs

```@docs
AbstractInterpolation
LinearInterpolation
QuadraticInterpolation
CubicSplineInterpolation
```

### Constructors

```@docs
LinearInterpolation(u::AbstractVector, t::AbstractVector)
QuadraticInterpolation(u::AbstractVector, t::AbstractVector)
CubicSplineInterpolation(u::AbstractVector, t::AbstractVector)
```

### Functions

```@docs
munge_data
_interp
```

## Einstein summation convention

```@docs
@einsum
```

## Curve fitting

### Structs

```@docs
OnceDifferentiable
OptimizationResults
LsqFitResult
```

### Constructors

```@docs
OnceDifferentiable(𝑓, p0::AbstractArray, 𝐹::AbstractArray)
```

### Functions

```@docs
value
value!
jacobian
jacobian!
levenberg_marquardt
curve_fit
```
