#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/04/27
#

#=
*Remarks* : Default Models

Flat model:

```math
\begin{equation}
m(\omega) = const.
\end{equation}
```

Gaussian model:

```math
\begin{equation}
m(\omega) = \frac{1}{\Gamma \sqrt{\pi}}
\exp\left[-\left(\frac{\omega}{\Gamma}\right)^2\right]
\end{equation}
```

Two Gaussians model:

```math
\begin{equation}
m(\omega) = 
\frac{1}{\Gamma \sqrt{\pi}}
\exp\left[-\left(\frac{\omega - \omega_0}{\Gamma}\right)^2\right]
+
\frac{1}{\Gamma \sqrt{\pi}}
\exp\left[-\left(\frac{\omega + \omega_0}{\Gamma}\right)^2\right]
\end{equation}
```

Lorentzian model:

```math
\begin{equation}
m(\omega)=\frac{\Gamma}{\pi(\Gamma^2+\omega^2)}
\end{equation}
```

Two Lorentzians model:

```math
\begin{equation}
m(\omega)=
\frac{\Gamma}{\pi[\Gamma^2+(\omega-\omega_0)^2]} +
\frac{\Gamma}{\pi[\Gamma^2+(\omega+\omega_0)^2]}
\end{equation}
```

Rise-And-Decay model:

```math
\begin{equation}
m(ω) = Γ\omega^2 \exp{(-\Gamma\omega)},\quad \omega \ge 0
\end{equation}
```
=#

"""
    build_flat_model(am::AbstractMesh)

Try to build a flat model in `am`. Then this model function is normalized.

See also: [`AbstractMesh`](@ref).
"""
function build_flat_model(am::AbstractMesh)
    model = ones(F64, length(am))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_gaussian_model(am::AbstractMesh, Γ::F64 = 2.0)

Try to build a Gaussian model, which is then normalized. The argument
`Γ` is used to control the width of the Gaussian peak.

See also: [`AbstractMesh`](@ref).
"""
function build_gaussian_model(am::AbstractMesh, Γ::F64 = 2.0)
    model = exp.(-(am.mesh / Γ) .^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_lorentzian_model(am::AbstractMesh, Γ::F64 = 2.0)

Try to build a Lorentzian model, which is then normalized. The argument
`Γ` is used to control the width of the Lorentzian peak.

See also: [`AbstractMesh`](@ref).
"""
function build_lorentzian_model(am::AbstractMesh, Γ::F64 = 2.0)
    model = (Γ / π) ./ ( Γ ^ 2.0 .+ (am.mesh) .^ 2.0 )
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_file_model(am::AbstractMesh, fn::String = "model.data")

Try to read a model function from external file (specified by `fn`). Note
that the mesh used to generate the model function must be compatible with
`am`. In addition, the model function will not be normalized.

See also: [`AbstractMesh`](@ref).
"""
function build_file_model(am::AbstractMesh, fn::String = "model.data")
    model = zeros(F64, length(am))
    open(fn, "r") do fin
        for i in eachindex(model)
            arr = parse.(F64, line_to_array(fin)[1:2])
            _mesh = arr[1]
            @assert abs(am[i] - _mesh) < 1e-6
            model[i] = arr[2]
        end
    end
    return model
end

"""
    build_func_model(fun::Function, am::AbstractMesh, kwargs...)

Try to build a model function by customized function `fun`. `kwargs`
denotes the arguments required by `fun`.

See also: [`AbstractMesh`](@ref).
"""
function build_func_model(fun::Function, am::AbstractMesh, kwargs...)
    model = fun.(am, kwargs...)
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end
