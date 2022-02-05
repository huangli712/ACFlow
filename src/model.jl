#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

function build_flat_model(am::AbstractMesh)
    model = ones(F64, length(am))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

#=
```math
\begin{equation}
m(\omega) = \frac{1}{\Gamma \sqrt{\pi}} \exp\left[-(\omega/\Gamma)^2\right]
\end{equation}
```
=#

function build_gaussian_model(am::AbstractMesh, Γ::F64 = 2.0)
    model = exp.(-(am.mesh / Γ) .^ 2.0) / (Γ * sqrt(π))
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end

function build_func_model(f::Function, am::AbstractMesh, kwargs...)
    model = f.(am, kwargs...)
    norm = dot(am.weight, model)
    model = model ./ norm
    return model
end
