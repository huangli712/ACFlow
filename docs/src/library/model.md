*For default model functions.*

The ACFlow toolkit supports various model functions, such as flat, Gaussian, Lorentzian, and a few unusual models. They are useful for the `MaxEnt` and `StochAC` solvers. In order to build these model functions, we need some additional parameters, including ``\Gamma``, ``s``, ``s_1``, and ``s_2``. They should be setup by using the parameter [`pmodel`](@ref pmodel).

```@index
Pages = ["model.md"]
```

## Flat Model

```@docs
build_flat_model
```

## Gaussian Models

!!! note

    This class includes the standard Gaussian model, shifted Gaussian model, and two-Gaussians model.

```@docs
build_gaussian_model
build_1gaussian_model
build_2gaussians_model
```

## Lorentzian Models

!!! note

    This class includes the standard Lorentzian model, shifted Lorentzian model, and two-Lorentzians model.

```@docs
build_lorentzian_model
build_1lorentzian_model
build_2lorentzians_model
```

## Unusual Models

!!! note

    This class includes the Rise-And-Decay model, file model, and function model.

```@docs
build_risedecay_model
build_file_model
build_func_model
```
