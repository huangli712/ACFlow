#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/08
#

#=
*Remarks* :

Flat model:

```math
\begin{equation}
m(\omega) = const.
\end{equation}
```

Gaussian model:

```math
\begin{equation}
m(\omega) = \frac{1}{\Gamma \sqrt{\pi}} \exp\left[-(\omega/\Gamma)^2\right]
\end{equation}
```
=#

"""
    build_flat_model(am::AbstractMesh)

Try to build a flat model, which is then normalized.

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

Try to build a gaussian model, which is then normalized. The parameter
`Γ` is used to control the width of the gaussian peak.

See also: [`AbstractMesh`](@ref).
"""
function build_gaussian_model(am::AbstractMesh, Γ::F64 = 2.0)
    model = exp.(-(am.mesh / Γ) .^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

"""
    build_file_model(am::AbstractMesh, fn::AbstractString)

Try to read a model function from external file (specified by `fn`). Note
that the mesh used to generate the model function must be compatible with
`am`. In addition, the model function will not be normalized.

See also: [`AbstractMesh`](@ref).
"""
function build_file_model(am::AbstractMesh, fn::AbstractString)
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
    build_func_model(f::Function, am::AbstractMesh, kwargs...)

Try to build a model function by customized function `f`. `kwargs` denotes
the parameters used in `f`.

See also: [`AbstractMesh`](@ref).
"""
function build_func_model(f::Function, am::AbstractMesh, kwargs...)
    model = f.(am, kwargs...)
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end
