*Define some essential mathematical functions.*

```@index
Pages = ["math.md"]
```

## Root Finding

```@docs
secant
newton
```

## Numerical Integrations

```@docs
trapz
simpson
```

## Numerical Derivatives

```@docs
second_derivative
gradient_via_fd
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

## Einstein Summation Convention

```@docs
@einsum
```

## Curve Fitting

### Structs

```@docs
OnceDifferentiable
LMOptimizationResults
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

## Numerical Optimization

### Structs

```@docs
BFGSDifferentiable
BFGSState
BFGSOptimizationResults
```

### Constructors

```@docs
BFGSDifferentiable(f, df, x::AbstractArray)
```

### Functions

```@docs
value(obj::BFGSDifferentiable)
gradient
value_gradient!
maxdiff
eval_δf
eval_Δf
eval_δx
eval_Δx
eval_resid
optimize
init_state
update_state!
update_g!
update_h!
trace
linesearch!
converged
```

## Line Search

### Structs

```@docs
LineSearchException
```

### Functions

```@docs
LS
```
