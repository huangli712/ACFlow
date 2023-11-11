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

```@docs
LinearInterpolation(u::AbstractVector, t::AbstractVector)
QuadraticInterpolation(u::AbstractVector, t::AbstractVector)
CubicSplineInterpolation(u::AbstractVector, t::AbstractVector)
```

## Einstein Summation Convention

```@docs
@einsum
```

## Curve Fitting

```@docs
curve_fit
```

## Numerical Optimization

```@docs
optimize
```