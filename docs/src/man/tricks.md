# [Tricks and tips](@id tricks)

*A comprehensive guide about how to perform reliable analytic continuation.*

As is well-known, analytic continuation is a tedious, tricky, and complicated job. This page contains all of the details (and experiences), that need to be known by the users, about how to get reliable and accurate spectra by using various solvers as implemented in the ACFlow toolkit.

!!! warning

    The following tricks and tips are **MY** personal viewpoints or experiences. They are not general rules. They might be wrong or ineffective or useless for some cases. They are definitely not suitable for the other analytic continuation codes.

## Contents

```@contents
Pages = ["tricks.md"]
Depth = 3
```

## [General rules](@id general)

1. It would be better to perform analytic continuation in Matsubara frequency axis, instead of imaginary time axis. See [`grid`](@ref grid) and [`ngrid`](@ref ngrid).

2. Employ the MaxEnt solver for a fast analytic continuation task. And then the stochastic methods can be used to get a better spectrum. See [`solver`](@ref solver).

3. If the input data is broken or discontinuous, please setup the [`grid`](@ref grid) parameter correctly.

## [MaxEnt solver](@id maxent)

1. The `chi2kink` and `bryan` algorithms are recommended. See [`method`](@ref maxent_method).

2. The Shannon-Jaynes entropy is recommented. But sometimes, if sharp features are essential, please choose the Bayesian Reconstruction entropy. See [`stype`](@ref stype).

## [NevanAC solver](@id nevanac)

1. It is extremely sensitive to the noise. So please make sure that the input data is smooth and is free of noise.

2. to be done

## [StochAC solver](@id stochac)

1. Increase `nfine` to 20000.

2. to be done.

## [StochSK solver](@id stochsk)

1. The `chi2min` algorithm is recommended. See [`method`](@ref stochsk_method).

2. to be done.

## [StochOM solver](@id stochom)

1. It is more efficient for Matsubara frequency Green's functions.

2. to be done.

## [StochPX solver](@id stochpx)

1. If the spectrum is expected to be broad, please set `method = 'mean'`. If the spectrum is expected to be ``\delta``-like, please set `method = 'best'`. See [`method`](@ref stochpx_method).

2. Run it parallelly (use `Pacrun.jl`).
