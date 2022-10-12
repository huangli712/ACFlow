# Models

*Define some type aliases and string constants for the ACFlow package.*

## Contents

```@contents
Pages = ["model.md"]
Depth = 2
```

## Index

```@index
Pages = ["model.md"]
```

## Flat model

```@docs
build_flat_model
```

## Gaussian models

!!! note

    The parameters ``\Gamma``, ``s``, ``s_1``, and ``s_2``, which are essential in defining the Gaussian-like models are setup by using the parameter [`pmesh`](@ref pmesh).

```@docs
build_gaussian_model
build_1gaussian_model
build_2gaussians_model
```

## Lorentzian models

```@docs
build_lorentzian_model
build_1lorentzian_model
build_2lorentzians_model
```

## Unusual models

```@docs
build_risedecay_model
build_file_model
build_func_model
```
