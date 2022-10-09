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

### Macros

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
