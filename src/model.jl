#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/08/26
#

#=
*Remarks* : *Default Models*

Now `ACFlow` supports various default model functions `m(ω)`. Note that
the `BarRat` solver (based on the barycentric rational function method),
the `NevanAC` solver (based on the Nevanlinna analytical continuation),
the `StochOM` solver (based on the stochastic optimization method), the
`StochSK` solver (based on the stochastic analytic continuation method),
and the `StochPX` solver (based on the stochastic pole expansion method)
do not need any default model functions. However, the `StochAC` solver
(based on the stochastic analytic continuation method as well) only
supports the `flat` default model function. The users can use the `mtype`
parameter to customize the default model that should be used.

These default model functions are summaried as follows.

* Flat model (keyword : `flat`):

```math
\begin{equation}
m(\omega) = const.
\end{equation}
```

* Gaussian model (keyword : `gauss`):

```math
\begin{equation}
m(\omega) = \frac{1}{\Gamma \sqrt{\pi}}
\exp\left[-\left(\frac{\omega}{\Gamma}\right)^2\right]
\end{equation}
```

* Shifted Gaussian model (keyword : `1gauss`):

```math
\begin{equation}
m(\omega) =
    \frac{1}{\Gamma \sqrt{\pi}}
    \exp\left[-\left(\frac{\omega - \omega_1}{\Gamma}\right)^2\right]
\end{equation}
```

* Two Gaussians model (keyword : `2gauss`):

```math
\begin{equation}
m(\omega) =
    \frac{1}{\Gamma \sqrt{\pi}}
    \exp\left[-\left(\frac{\omega - \omega_1}{\Gamma}\right)^2\right]
    +
    \frac{1}{\Gamma \sqrt{\pi}}
    \exp\left[-\left(\frac{\omega - \omega_2}{\Gamma}\right)^2\right]
\end{equation}
```

* Lorentzian model (keyword : `lorentz`):

```math
\begin{equation}
m(\omega)=\frac{\Gamma}{\pi(\Gamma^2+\omega^2)}
\end{equation}
```

* Shifted Lorentzian model (keyword : `1lorentz`):

```math
\begin{equation}
m(\omega) =
    \frac{\Gamma}{\pi[\Gamma^2+(\omega - \omega_1)^2]}
\end{equation}
```

* Two Lorentzians model (keyword : `2lorentz`):

```math
\begin{equation}
m(\omega) =
    \frac{\Gamma}{\pi[\Gamma^2+(\omega - \omega_1)^2]} +
    \frac{\Gamma}{\pi[\Gamma^2+(\omega - \omega_2)^2]}
\end{equation}
```

* Rise-And-Decay model (keyword : `risedecay`):

```math
\begin{equation}
m(ω) = Γ\omega^2 \exp{(-\Gamma\omega)},\quad \omega \ge 0
\end{equation}
```

Be careful, ``\Gamma``, ``\omega_1``, and ``\omega_2`` are the model
parameters. Their default values are defined in `make_model()`, but
you can modified them according to the `pmodel` parameter.

---
=#

"""
    build_flat_model(am::AbstractMesh)

Try to build a flat model in `am`. Then this model function is normalized.
Only this model function is suitable for the `StochAC` solver.

### Arguments
* am -> Real frequency mesh.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_flat_model(am::AbstractMesh)
    model = ones(F64, length(am))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_gaussian_model(am::AbstractMesh, Γ::F64)

Try to build a Gaussian model, which is then normalized. The argument
`Γ` is used to control the width of the Gaussian peak.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control width of the peak.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_gaussian_model(am::AbstractMesh, Γ::F64)
    model = exp.(-(am.mesh / Γ) .^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_1gaussian_model(am::AbstractMesh, Γ::F64, s::F64)

Try to build a shifted Gaussian model, which is then normalized. The
argument `Γ` is used to control the width of the Gaussian peak, and `s`
means the shift of the central peak. If `s > 0`, the peak is shifted to
positive half-axis, and vice versa.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control width of the gaussian peak.
* s -> Distance of the peak from ω = 0.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_1gaussian_model(am::AbstractMesh, Γ::F64, s::F64)
    model = similar(am.mesh)
    @. model = exp(-((am.mesh - s) / Γ) ^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_2gaussians_model(am::AbstractMesh, Γ::F64, s₁::F64, s₂::F64)

Try to build a Two Gaussians model, which is then normalized. The
argument `Γ` is used to control the width of the Gaussian peak, and
`s₁` and `s₂` denote the centers of the two peaks.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control width of the gaussian peaks.
* s₁ -> Distance of the first peak from ω = 0.
* s₂ -> Distance of the second peak from ω = 0.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_2gaussians_model(am::AbstractMesh, Γ::F64, s₁::F64, s₂::F64)
    model = similar(am.mesh)
    @. model = exp(-((am.mesh - s₁) / Γ) ^ 2.0) / (Γ * sqrt(π))
    @. model += exp(-((am.mesh - s₂) / Γ) ^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_lorentzian_model(am::AbstractMesh, Γ::F64)

Try to build a Lorentzian model, which is then normalized. The argument
`Γ` is used to control the width of the Lorentzian peak.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control broadening of the lorentzian peak.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_lorentzian_model(am::AbstractMesh, Γ::F64)
    model = (Γ / π) ./ ( Γ ^ 2.0 .+ (am.mesh) .^ 2.0 )
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_1lorentzian_model(am::AbstractMesh, Γ::F64, s::F64)

Try to build a shifted Lorentzian model, which is then normalized. The
argument `Γ` is used to control the width of the Lorentzian peak, and `s`
means the shift of the central peak. If `s > 0`, the peak is shifted to
positive half-axis, and vice versa.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control broadening of the lorentzian peak.
* s -> Distance of the peak from ω = 0.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_1lorentzian_model(am::AbstractMesh, Γ::F64, s::F64)
    model = (Γ / π) ./ ( Γ ^ 2.0 .+ (am.mesh - s) .^ 2.0 )
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_2lorentzians_model(am::AbstractMesh, Γ::F64, s₁::F64, s₂::F64)

Try to build a Two-Lorentzians model, which is then normalized. The
argument `Γ` is used to control the width of the Lorentzian peak, and
`s₁` and `s₂` denote the centers of the two peaks.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control broadening of the lorentzian peaks.
* s₁ -> Distance of the first peak from ω = 0.
* s₂ -> Distance of the second peak from ω = 0.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_2lorentzians_model(am::AbstractMesh, Γ::F64, s₁::F64, s₂::F64)
    model = similar(am.mesh)
    @. model = (Γ / π) / ( Γ ^ 2.0 + (am.mesh - s₁) ^ 2.0 )
    @. model += (Γ / π) / ( Γ ^ 2.0 + (am.mesh - s₂) ^ 2.0 )
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_risedecay_model(am::AbstractMesh, Γ::F64)

Try to build a Rise-And-Decay model, which is then normalized. This model
function is defined on positive half-axis, so it is more suitable for the
bosonic response function.

### Arguments
* am -> Real frequency mesh.
* Γ -> Parameter to control rise and decay of the peak.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_risedecay_model(am::AbstractMesh, Γ::F64)
    @assert am[1] ≥ 0.0
    model = Γ * (am.mesh .^ 2.0) .* exp.(-Γ * am.mesh)
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_file_model(am::AbstractMesh, fn::String)

Try to read a model function from external file (specified by `fn`). Note
that the mesh used to generate the model function must be compatible with
`am`. In addition, the model function will be normalized automatically.

### Arguments
* am -> Real frequency mesh.
* fn -> Filename for external model function.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_file_model(am::AbstractMesh, fn::String)
    # Allocate memory
    model = zeros(F64, length(am))
    #
    # Read model function
    open(fn, "r") do fin
        for i in eachindex(model)
            arr = parse.(F64, line_to_array(fin)[1:2])
            _mesh = arr[1]
            @assert abs(am[i] - _mesh) < 1e-6
            model[i] = arr[2]
        end
    end
    #
    # Normalize model function
    norm = dot(am.weight, model)
    model = model ./ norm

    return model
end

"""
    build_func_model(fun::Function, am::AbstractMesh, kwargs...)

Try to build a model function by customized function `fun`. `kwargs`
denotes the arguments required by `fun`. Actually, this feature does
**not** really work.

### Arguments
* fun -> A external funtion call to build the model function.
* am -> Real frequency mesh.
* kwargs -> Arguments that passed to `fun`.

### Returns
* model -> Default model function, m(ω). ω is compatible with `am`.

See also: [`AbstractMesh`](@ref).
"""
function build_func_model(fun::Function, am::AbstractMesh, kwargs...)
    model = fun.(am, kwargs...)
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end
