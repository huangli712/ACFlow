using LinearAlgebra
using DoubleFloats
using DelimitedFiles
using Printf
using Plots

"""
    OperatorType Bose Fermi

Operator type of physical observables
"""
@enum OperatorType Bose Fermi

"""
    Identity matrix
"""
function eye(dtype::DataType, n::Int64)
    return Matrix{dtype}(I,n,n)
end

"""
    mobius_transform(z, Y)

Mobius_transform transforms the upper half complex plane to the unit disk
in the complex plane, in which point `Y` is mapped to the center of `D`.
"""
function mobius_transform(z, Y) 
    if abs(Y) == 0 @warn "Y=0.0 maps all points to 0.0!" end
    if imag(z) < 0 @warn "Im(z) ≥ 0 is required."  end
    if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required."  end
    num = z - Y
    den = z - Y'
    return num/den
end

"""
    inverse_mobius_transform(z, Y)
    
Inverse Mobius transform transforms the unit disk in the complex plane
to the upper half complex plane.
"""
function inverse_mobius_transform(z, Y)
    if abs(z) ≥ 1 @warn "|z|<1 is required." end
    if abs(Y) == 0 @warn "Y=0.0 maps all points to 0!" end
    if imag(Y) < 0 @warn "Im(Y) ≥ 0 is required." end
    dtype = eltype(z)
    num = Y - z*Y'
    den = one(dtype) - z
    return num/den
end

"""
    pick_matrix(x, y)

Pick matrix of initial data `{z, f(z)}` with in a unit cell.
"""
function pick_matrix(x::AbstractVector, y::AbstractVector)
    if (all(abs.(x) .> 0) && all(abs.(y) .> 0)) == false
        @error "pick_matrix: initial data and target data should be in a unit cell"
    elseif length(x) != length(y)
        @error DimensionMismatch
    else
        N = length(y)
        pick = similar(y, N, N)
        for i = 1: N, j = 1: N
            num = 1 - y[i] * y[j]'
            den = 1 - x[i] * x[j]'
            pick[i,j] = num/den
        end
        return pick
    end
end

"""
    isNevanlinnasolvable(x, y[;tolerance=1.e-10])
    
Check if the initial data `{z, f(z)}` can be interpolated by generalized
schur algorithm for Nevanlinna functions.
"""
function isNevanlinnasolvable(x::AbstractVector, y::AbstractVector; 
    tolerance::AbstractFloat = 1.e-10)
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @error "Target data should be in the upper half complex plane"
    else
        x = mobius_transform.(x, oneunit(x[1])im)
        y = mobius_transform.(y, oneunit(y[1])im)
        pick = pick_matrix(x, y)
        evals = eigvals(pick) 
        return all((real.(evals) .+ tolerance) .>= 0.), minimum(real.(evals))
    end
end

#=
```math
θ(n-1) = [A[1,1]*θ(n) + A[1,2]]/[A[2,1]*θ(n)+A[2,2]]
```
=#

"""
    _recursion(A::AbstractMatrix, θn::Number) 

Recursion relation of contractive functions `θ(n)` and `θ(n-1)`,the relation is:
"""
function _recursion(A::AbstractMatrix, θn::Number)
    num = A[1,1] * θn + A[1,2]
    den = A[2,1] * θn + A[2,2]
    return num/den
end

#=
```math
θ(n+1) = [-A[2,2]*θ(n) + A[1,2]]/[A[2,1]*θ(n)-A[1,1]]
```
=#

"""
    _inv_recursion(A::AbstractMatrix, θn::Number) 

Recursion relation of contractive functions `θ(n)` and `θ(n+1)`,the relation is:
'''
"""
function _inv_recursion(A::AbstractMatrix, θp::Number)
    num = -A[2,2] * θp + A[1,2]
    den = A[2,1] * θp - A[1,1]
    return num/den
end

"""
    _coefficient(z, xj, ϕj) 

Calculate the coefficient of the recursion relation in generalized_schur algorithm.
"""
function coefficient(z::Number, xj::Number, ϕj::Number) 
    A = zeros(typeof(z),2,2)
    A[1,1] = mobius_transform(z, xj)
    A[1,2] = ϕj
    A[2,1] = ϕj' * mobius_transform(z, xj)
    A[2,2] = one(typeof(z))
    return A
end

"""
    schur_parameter([ftype::DataType = T,] x::AbstractVector{T}, y::AbstractVector{T}) where T
    
Evaluate Schur parameters for contractive functions `y(x)` within a unit
circle, return a list of Schur parameters
"""
function schur_parameter(x::AbstractVector{T}, y::AbstractVector{T}) where T
    M = length(y)
    ϕ = zeros(T, M); ϕ[1] = y[1]
    abcd = fill(eye(T,2), M)
    factor = fill(zeros(T,2,2), M-1)
    for j = 1:(M-1)
        for k=j:M
            prod = coefficient(x[k], x[j], ϕ[j])
            abcd[k] *= prod
        end
        ϕ[j+1] = _inv_recursion(abcd[j+1], y[j+1])
        factor[j] = abcd[j+1]
    end
    return ϕ
end

"""
    generalized_schur(z::Number,
                      x::AbstractVector{T},
                      y::AbstractVector{T}[; init_func::Function = x->zero(T)]) where T

The generalized Schur algorithm that extrapolates beween `{x,y}` and
generate a contractive function `f(z)`, return its value at `z`.
"""
function generalized_schur(z::Number,
                           x::AbstractVector,
                           y::AbstractVector;
                           init_func::Function = z -> zero(eltype(y))) 
    if all(imag.(x) .≥ 0) == false 
        @error "Initial data should be in the upper half complex plane"
    elseif all(abs.(y) .≤ 1) == false 
        @error "Target data should be in the unit circle"
    else
        M = length(y)
        ϕ = schur_parameter(x,y)
        abcd = eye(eltype(y),2)
        for j = 1:M
            abcd *= coefficient(z,x[j],ϕ[j])
        end
        return _recursion(abcd, init_func(z))
    end
end

"""
    nevanlinna(z::Number,
               x::AbstractVector{T},
               y::AbstractVector{T}[; init_func::Function = x->zero(T)]) where T

The Nevanlinna Interpolation algorithm that extrapolates beween `{x,y}`
and generate a nevanlinna function `f(z)`, return its value at `z`.    
"""
function nevanlinna(z::Number,
                    x::AbstractVector{T},
                    y::AbstractVector{T}; 
                    init_func::Function = z -> zero(T)) where T
    if all(imag.(x) .≥ 0) == false 
        @warn "Initial data should be in the upper half complex plane"
    elseif all(imag.(y) .≥ 0) == false
        @warn "Target data should be in the upper half complex plane"
    end
    y = mobius_transform.(y, oneunit(y[1])im)
    res = generalized_schur(z, x, y; init_func)
    return inverse_mobius_transform(res, oneunit(res)im)
end

"""
    toNevanlinnadata(x::AbstractVector,
                     y::AbstractVector,
                     type::OperatorType[, float_type::DataType = Double64])

Convert Masubara frequency correlation function data to Nevanlinna data
"""
function toNevanlinnadata(x::AbstractVector,
                          y::AbstractVector,
                          operator_type::OperatorType,
                          float_type::DataType = Double64)
    complex_type = Complex{float_type}
    Im = one(float_type)im
    x = convert(Vector{complex_type}, x)
    y = convert(Vector{complex_type}, y)
    x = Im * x
    if operator_type == Fermi
        y = -y
    else
        y = -x .* y
    end
    return x, y
end

"""
    spectral_function(operator_type::OperatorType,
                        x::AbstractVector,
                        y::AbstractVector[, float_type::DataType=Double64];
                        η::Real = 1.e-5,
                        init_func::Function= z -> 0.0)

Calculate the spectral function `A(ω)` for given dataset `{x=ωn, y=G(iωn)}` at `ω`.
"""
function spectral_function(operator_type::OperatorType,
                            x::AbstractVector, y::AbstractVector,
                            float_type::DataType=Double64;
                            η::Real = 1.e-5,
                            init_func::Function= z -> zero(Complex{float_type}))
    x, y = toNevanlinnadata(x, y, operator_type, float_type)
    if isNevanlinnasolvable(x,y)[1] == false @warn "Nevanlinna unsolvable!" end

    ωmax = 4π
    Nω = 500
    L = 2*ωmax
    ωlist = convert(Vector{float_type}, (-Nω/2:Nω/2-1)*L/Nω)

    Alist = []
    for i in eachindex(ωlist) 
        ω = convert(float_type, ωlist[i])
        z = ω + one(float_type)im * η 
        res = 2*imag(nevanlinna(z, x, y; init_func))
        
        if operator_type == Fermi
            push!(Alist, res)
        else
            push!(Alist, res/ω)
        end
    end

    return ωlist, Alist
end

"""
    Read Masubara Green's function data
"""
function readGF(float_type::DataType, path::String; rev::Bool=true, num::Integer = 100)
    Im = one(float_type)im  #define imaginary unit
    ctype =  Complex{float_type} 
    d = readdlm(path)
    x = ctype.(d[:,1]) 
    y = float_type.(d[:,2]) + Im * float_type.(d[:,3])
    num = min(length(x), num)
    x1, y1 = x[1:num], y[1:num]
    if rev == true
        x1 = reverse(x1); y1=reverse(y1)
    end
    return x1, y1
end

#A1data = readdlm("./data/A_delta_eta_0.05.txt")
#ωlist = A1data[:,1]
#A1 = A1data[:, 2]
x1, y1 = readGF(Float64, "repr.data")
wlist, A1ac = spectral_function(Fermi, x1, y1)
p1 = plot(wlist, A1ac, line=2, label="exact")
#plot!(wlist, A1ac, line=(1,:dash), marker=2, label="ac")
