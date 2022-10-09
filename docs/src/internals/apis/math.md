# Math

*Define some essential mathematical functions.*

## Root finding

### Functions

```@docs
secant
newton
```

## Numerical integrations

### Functions

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

```@docs
LinearInterpolation(u::AbstractVector, t::AbstractVector)
```


## Functions

```@docs
@einsum
curve_fit
levenberg_marquardt
value
value!
jacobian
jacobian!
munge_data
_interp
```

```@docs
LsqFitResult
OptimizationResults
OnceDifferentiable
```