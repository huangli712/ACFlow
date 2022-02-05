include("global.jl")
include("mesh.jl")
include("util.jl")
include("model.jl")

using Distributed
using Dates
using LinearAlgebra

wmin = -5.0
wmax = 5.0
nmesh = 101

lm = LinearMesh(nmesh, wmin, wmax)
tm = TangentMesh(nmesh, wmin, wmax)

f(x) = x ^ 2.0
I(x) = x ^ 3.0 / 3.0

y_lm = f.(lm)
y_tm = f.(tm)

area_std = I(wmax) - I(wmin)

area_lm = area(lm, y_lm)
area_lm_trapz = trapz(lm.mesh, y_lm)
area_lm_trapz_linear = trapz(lm.mesh, y_lm, true)
area_lm_simpson = simpson(lm.mesh, y_lm)
area_tm = area(tm, y_tm)
area_tm_trapz = trapz(tm.mesh, y_tm)
#@show area_lm, area_lm_trapz, area_lm_trapz_linear, area_lm_simpson
#@show area_tm, area_tm_trapz
#@show area_std

Gmodel1 = build_gaussian_model(lm)

f(x, Γ1, Γ2) = exp(-(x / Γ1) ^ 2.0) / (Γ2 * sqrt(π))

Gmodel2 = build_func_model(f, lm, 3.0, 10.0)
#@show Gmodel1
#@show Gmodel2

#welcome()
#overview()

function myfun(x, b, c)
    return x ^ 2.0 + b * x + c
end

function myfun2(x, b)
    return x ^ 3.0 - b
end

s = secant(myfun2, 1.0, 8.0, maxiter = 40, tol = 1e-8)
println(s)