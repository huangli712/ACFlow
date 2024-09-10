#
# Project : Gardenia
# Source  : math.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/10
#

#
# BE CAREFUL:
#
# The algorithms as implemented here are designed for the ACFlow package
# only. Some of them are oversimplified to get better performance. Thus
# they won't work for general cases sometimes. Anyway, it is not a good
# idea to adopt these functions outside the ACFlow package. Please use
# popular Julia's packages in your codes.
#

#=
### *Math* : *Root Finding*
=#

"""
    secant(func, x0, args...)

It implements the well-known secant algorithm to locate root of a given
polynomial function. Here, `func` means the function, `x0` is the initial
guess, and `args...` denotes the required arguments for function call to
`func`. In addition, the maximum iterations and tolerance are controled
by `maxiter` and `tol`, respectively. Be careful, `func` must be a single
variable function.

### Arguments
See above explanations.

### Returns
* 𝑝 -> The solution.

See also: [`newton`](@ref).
"""
function secant(func, x0, args...; maxiter::I64 = 50, tol::F64 = 1.48e-8)
    eps = 1.0e-4
    funcalls = 0

    p0 = 1.0 * x0
    p1 = x0 * (1.0 + eps)
    if p1 ≥ 0.0
        p1 = p1 + eps
    else
        p1 = p1 - eps
    end

    q0 = func(p0, args...)
    funcalls = funcalls + 1
    q1 = func(p1, args...)
    funcalls = funcalls + 1

    if abs(q1) < abs(q0)
        p0, p1 = p1, p0
        q0, q1 = q1, q0
    end

    for _ = 1:maxiter
        if q1 == q0
            if p1 != p0
                error("Tolerance is reached in secant()!")
            end
            p = (p1 + p0) / 2.0
            return p
        else
            if abs(q1) > abs(q0)
                p = (-q0 / q1 * p1 + p0) / (1 - q0 / q1)
            else
                p = (-q1 / q0 * p0 + p1) / (1 - q1 / q0)
            end
        end

        if abs(p - p1) < tol
            return p
        end

        p0, q0 = p1, q1
        p1 = p
        q1 = func(p1, args...)
        funcalls = funcalls + 1
    end
end

"""
    newton(
        fun::Function,
        guess,
        kwargs...;
        maxiter::I64 = 20000,
        mixing::F64 = 0.5
    )

It implements the well-known newton algorithm to locate root of a given
polynomial function. Here, `fun` means the function, `guess` is the initial
solution, and `kwargs...` denotes the required arguments for `fun`. Please
be careful, `func` is a multiple variable function. It not only returns
the value, but also the jacobian matrix of the function.

### Arguments
See above explanations.

### Returns
* sol -> Solution.
* call -> Counter for function call to `fun()`.

See also: [`secant`](@ref).
"""
function newton(
    fun::Function,
    guess,
    kwargs...;
    maxiter::I64 = 20000,
    mixing::F64 = 0.5
    )
    function _apply(feed::Vector{T}, f::Vector{T}, J::Matrix{T}) where {T}
        resid = nothing
        step = 1.0
        limit = 1e-4
        try
            resid = - pinv(J) * f
        catch
            resid = zeros(F64, length(feed))
        end
        if any(x -> x > limit, abs.(feed))
            ratio = abs.(resid ./ feed)
            max_ratio = maximum( ratio[ abs.(feed) .> limit ] )
            if max_ratio > 1.0
                step = 1.0 / max_ratio
            end
        end
        return feed + step .* resid
    end

    counter = 0
    feeds = []
    backs = []

    f, J = fun(guess, kwargs...)
    back = _apply(guess, f, J)
    push!(feeds, guess)
    push!(backs, back)

    while true
        counter = counter + 1
        feed = feeds[end] + mixing * (backs[end] - feeds[end])

        f, J = fun(feed, kwargs...)
        back = _apply(feed, f, J)
        push!(feeds, feed)
        push!(backs, back)

        any(isnan.(back)) && error("Got NaN!")
        if counter > maxiter || maximum( abs.(back - feed) ) < 1.e-4
            break
        end
    end

    counter > maxiter && error("Tolerance is reached in newton()!")

    return back, counter
end

#=
### *Math* : *Numerical Integrations*
=#

"""
    trapz(x::AbstractMesh, y::AbstractVector{T}) where {T<:N64}

Perform numerical integration by using the composite trapezoidal rule.

### Arguments
* x -> Real frequency mesh.
* y -> Function values at real axis.

### Returns
* ℐ -> The final value.

See also: [`simpson`](@ref).
"""
function trapz(x::AbstractMesh, y::AbstractVector{T}) where {T<:N64}
    value = dot(x.weight, y)
    return value
end

"""
    trapz(
        x::AbstractVector{S},
        y::AbstractVector{T},
        linear::Bool = false
    ) where {S<:Number, T<:Number}

Perform numerical integration by using the composite trapezoidal rule.
Note that it supports arbitrary precision via BigFloat.

### Arguments
* x      -> Real frequency mesh.
* y      -> Function values at real axis.
* linear -> Whether the given mesh is linear?

### Returns
* ℐ -> The final value.

See also: [`simpson`](@ref).
"""
function trapz(
    x::AbstractVector{S},
    y::AbstractVector{T},
    linear::Bool = false
    ) where {S<:Number, T<:Number}
    # For linear mesh
    if linear
        h = x[2] - x[1]
        value = y[1] + y[end] + 2.0 * sum(y[2:end-1])
        value = h * value / 2.0
    # For non-equidistant mesh
    else
        len = length(x)
        dx = view(x, 2:len) .- view(x, 1:(len-1))
        y_forward = view(y, 2:len)
        y_backward = view(y, 1:(len-1))
        value = sum(0.5 * (y_forward .+ y_backward) .* dx)
    end

    return value
end

"""
    simpson(
        x::AbstractVector{S},
        y::AbstractVector{T}
    ) where {S<:Number, T<:Number}

Perform numerical integration by using the simpson rule. Note that the
length of `x` and `y` must be odd numbers. And `x` must be a linear and
uniform mesh.

### Arguments
* x -> Real frequency mesh.
* y -> Function values at real axis.

### Returns
* ℐ -> The final value.

See also: [`trapz`](@ref).
"""
function simpson(
    x::AbstractVector{S},
    y::AbstractVector{T}
    ) where {S<:Number, T<:Number}
    h = (x[2] - x[1]) / 3.0

    even_sum = 0.0
    odd_sum = 0.0
    for i = 2:length(x)-1
        if iseven(i)
            even_sum = even_sum + y[i]
        else
            odd_sum = odd_sum + y[i]
        end
    end

    return h * (y[1] + y[end] + 4.0 * even_sum + 2.0 * odd_sum)
end

#=
### *Math* : *Numerical Derivatives*
=#

"""
    second_derivative(x::AbstractVector, y::AbstractVector)

Compute second derivative y''(x). If the length of `x` and `y` is `N`,
the length of the returned vector is `N-2`.

### Arguments
* x -> Real frequency mesh.
* y -> Function values at real axis.

### Returns
* val -> y''(x).
"""
function second_derivative(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y)

    N = length(x)
    h₁ = view(x, 2:(N-1)) - view(x, 1:(N-2))
    h₂ = view(x, 3:N) - view(x, 2:(N-1))

    y_forward = view(y, 3:N)
    y_mid = view(y, 2:(N-1))
    y_backward = view(y, 1:(N-2))

    num = h₁ .* y_forward + h₂ .* y_backward - (h₁ + h₂) .* y_mid
    den = (h₂.^2) .* h₁ + (h₁.^2) .* h₂
    return 2 .* num ./ den
end

"""
    gradient_via_fd(f, x)

Compute ∂f/∂x via finite difference method. It is less accurate and much
slower than the automatic differentiation approach. Actually, we won't
use this function to calculate gradient. The `Zygote.gradient()` function
is always a better choice.

### Arguments
* x -> Real frequency mesh.
* f -> Function values at real axis, f(x).

### Returns
* val -> ∂f/∂x.
"""
function gradient_via_fd(f, x)
    𝑠 = cbrt(eps(F64))

    ∇𝑓 = zero(x)
    𝑥 = copy(x)

    @inbounds for i in eachindex(x)
        ϵ = max(𝑠*abs(x[i]), 𝑠)
        #
        𝑥ᵢ = x[i]
        #
        𝑥[i] = 𝑥ᵢ + ϵ
        δ𝑓 = f(𝑥)
        𝑥[i] = 𝑥ᵢ - ϵ
        δ𝑓 -= f(𝑥)
        𝑥[i] = 𝑥ᵢ
        ∇𝑓[i] = real(δ𝑓 / (2 * ϵ))
        #
        𝑥[i] = 𝑥ᵢ + im * ϵ
        δ𝑓 = f(𝑥)
        𝑥[i] = 𝑥ᵢ - im * ϵ
        δ𝑓 -= f(𝑥)
        𝑥[i] = 𝑥ᵢ
        ∇𝑓[i] -= im * imag(δ𝑓 / (2 * im * ϵ))
    end
    ∇𝑓
end

#=
### *Math* : *Interpolations*
=#

#=
*Remarks* :

The following codes are used to perform interpolations. Three algorithms
are implemented. They are linear interpolation, quadratic interpolation,
and cubic spline interpolation. Note that these implementations are taken
directly from

* https://github.com/PumasAI/DataInterpolations.jl

Of cource, small modifications and simplifications are made.
=#

"""
    AbstractInterpolation

It represents an abstract interpolation engine, which is used to build
the internal type system.
"""
abstract type AbstractInterpolation{FT,T} <: AbstractVector{T} end

"""
    LinearInterpolation

It represents the linear interpolation algorithm.
"""
struct LinearInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    t::tType
    function LinearInterpolation{FT}(u,t) where {FT}
        new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
    end
end

"""
    LinearInterpolation(u::AbstractVector, t::AbstractVector)

Create a LinearInterpolation struct. Note that `u` and `t` denote `f(x)`
and `x`, respectively.
"""
function LinearInterpolation(u::AbstractVector, t::AbstractVector)
    u, t = munge_data(u, t)
    LinearInterpolation{true}(u,t)
end

"""
    QuadraticInterpolation

It represents the quadratic interpolation algorithm.
"""
struct QuadraticInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    t::tType
    function QuadraticInterpolation{FT}(u,t) where {FT}
        new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
    end
end

"""
    QuadraticInterpolation(u::AbstractVector, t::AbstractVector)

Create a QuadraticInterpolation struct. Note that `u` and `t` denote
`f(x)` and `x`, respectively.
"""
function QuadraticInterpolation(u::AbstractVector, t::AbstractVector)
    u, t = munge_data(u, t)
    QuadraticInterpolation{true}(u,t)
end

"""
    CubicSplineInterpolation

It represents the cubic spline interpolation algorithm.
"""
struct CubicSplineInterpolation{uType,tType,hType,zType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    t::tType
    h::hType
    z::zType
    function CubicSplineInterpolation{FT}(u,t,h,z) where {FT}
        new{typeof(u),typeof(t),typeof(h),typeof(z),FT,eltype(u)}(u,t,h,z)
    end
end

"""
    CubicSplineInterpolation(u::AbstractVector, t::AbstractVector)

Create a CubicSplineInterpolation struct. Note that `u` and `t` denote
`f(x)` and `x`, respectively.
"""
function CubicSplineInterpolation(u::AbstractVector, t::AbstractVector)
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k+1] - t[k], 1:length(t)-1), 0)
    dl = h[2:n+1]
    d_tmp = 2 .* (h[1:n+1] .+ h[2:n+2])
    du = h[2:n+1]
    tA = Tridiagonal(dl,d_tmp,du)
    d = map(i -> i == 1 || i == n + 1 ? 0 : 6(u[i+1] - u[i]) / h[i+1] - 6(u[i] - u[i-1]) / h[i], 1:n+1)
    z = tA\d
    CubicSplineInterpolation{true}(u,t,h[1:n+1],z)
end

"""
    munge_data(u::AbstractVector{<:Real}, t::AbstractVector{<:Real})

Preprocess the input data. Note that `u` and `t` denote `f(x)` and `x`,
respectively.
"""
function munge_data(u::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    return u, t
end

"""
    munge_data(u::AbstractVector, t::AbstractVector)

Preprocess the input data. Note that `u` and `t` denote `f(x)` and `x`,
respectively.
"""
function munge_data(u::AbstractVector, t::AbstractVector)
    Tu = Base.nonmissingtype(eltype(u))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == length(u)
    non_missing_indices = collect(
        i for i in 1:length(t) if !ismissing(u[i]) && !ismissing(t[i])
    )
    newu = Tu.([u[i] for i in non_missing_indices])
    newt = Tt.([t[i] for i in non_missing_indices])
    return newu, newt
end

"""
    (interp::AbstractInterpolation)(t::Number)

Support `interp(t)`-like function call, where `interp` is the instance of
any AbstractInterpolation struct and `t` means the point that we want to
get the interpolated value.
"""
(interp::AbstractInterpolation)(t::Number) = _interp(interp, t)

"""
    _interp(A::LinearInterpolation{<:AbstractVector}, t::Number)

To implement the linear interpolation algorithm.

See also: [`LinearInterpolation`](@ref).
"""
function _interp(A::LinearInterpolation{<:AbstractVector}, t::Number)
    idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
    θ = (t - A.t[idx])/(A.t[idx + 1] - A.t[idx])
    return (1 - θ)*A.u[idx] + θ*A.u[idx+1]
end

"""
    _interp(A::QuadraticInterpolation{<:AbstractVector}, t::Number)

To implement the quadratic interpolation algorithm.

See also: [`QuadraticInterpolation`](@ref).
"""
function _interp(A::QuadraticInterpolation{<:AbstractVector}, t::Number)
    idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 2))
    i₀, i₁, i₂ = idx, idx + 1, idx + 2
    l₀ = ((t - A.t[i₁])*(t - A.t[i₂]))/((A.t[i₀] - A.t[i₁])*(A.t[i₀] - A.t[i₂]))
    l₁ = ((t - A.t[i₀])*(t - A.t[i₂]))/((A.t[i₁] - A.t[i₀])*(A.t[i₁] - A.t[i₂]))
    l₂ = ((t - A.t[i₀])*(t - A.t[i₁]))/((A.t[i₂] - A.t[i₀])*(A.t[i₂] - A.t[i₁]))
    return A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

"""
    _interp(A::CubicSplineInterpolation{<:AbstractVector{<:Number}}, t::Number)

To implement the cubic spline interpolation algorithm.

See also: [`CubicSplineInterpolation`](@ref).
"""
function _interp(A::CubicSplineInterpolation{<:AbstractVector{<:Number}}, t::Number)
    i = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
    I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
    C = (A.u[i+1]/A.h[i+1] - A.z[i+1]*A.h[i+1]/6)*(t - A.t[i])
    D = (A.u[i]/A.h[i+1] - A.z[i]*A.h[i+1]/6)*(A.t[i+1] - t)
    I + C + D
end

export munge_data
export _interp

#=
### *Math* : *Einstein Summation Convention*
=#

#=
*Remarks* :

The following codes realize a single macro `@einsum`, which implements
*similar* notation to the Einstein summation convention to flexibly
specify operations on Julia `Array`s, similar to numpy's `einsum`
function (but more flexible!).

This implementation was borrowed from the following github repository:

* https://github.com/ahwillia/Einsum.jl

Of course, small modifications are made.
=#

"""
    @einsum(ex)

Perform Einstein summation like operations on Julia `Array`s.

### Examples

Basic matrix multiplication can be implemented as:

```julia
@einsum A[i, j] := B[i, k] * C[k, j]
```

If the destination array is preallocated, then use `=`:

```julia
A = ones(5, 6, 7) # Preallocated space
X = randn(5, 2)
Y = randn(6, 2)
Z = randn(7, 2)

# Store the result in A, overwriting as necessary:
@einsum A[i, j, k] = X[i, r] * Y[j, r] * Z[k, r]
```

If destination is not preallocated, then use `:=` to automatically create
a new array for the result:

```julia
X = randn(5, 2)
Y = randn(6, 2)
Z = randn(7, 2)

# Create new array B with appropriate dimensions:
@einsum B[i, j, k] := X[i, r] * Y[j, r] * Z[k, r]
```
"""
macro einsum(ex)
    _einsum(ex)
end

# Core function
function _einsum(expr::Expr, inbounds = true)
    # Get left hand side (lhs) and right hand side (rhs) of equation
    lhs = expr.args[1]
    rhs = expr.args[2]

    # Get info on the left hand side
    lhs_arrays, lhs_indices, lhs_axis_exprs = extractindices(lhs)

    # Get info on the right hand side
    rhs_arrays, rhs_indices, rhs_axis_exprs = extractindices(rhs)

    check_index_occurrence(lhs_indices, rhs_indices)

    # Remove duplicate indices on left hand and right hand side
    # and ensure that the array sizes match along these dimensions.
    dimension_checks = Expr[]

    # Remove duplicate indices on the right hand side
    for i in reverse(eachindex(rhs_indices))
        duplicated = false
        di = rhs_axis_exprs[i]

        for j = 1:(i - 1)
            if rhs_indices[j] == rhs_indices[i]
                duplicated = true
                dj = rhs_axis_exprs[j]
                push!(dimension_checks, :(@assert $dj == $di))
            end
        end

        for j = eachindex(lhs_indices)
            if lhs_indices[j] == rhs_indices[i]
                dj = lhs_axis_exprs[j]
                if Meta.isexpr(expr, :(:=))
                    # expr.head is :=
                    # Infer the size of the lhs array
                    lhs_axis_exprs[j] = di
                else
                    # expr.head is =, +=, *=, etc.
                    lhs_axis_exprs[j] = :(min($dj, $di))
                end
                duplicated = true
            end
        end

        if duplicated
            deleteat!(rhs_indices, i)
            deleteat!(rhs_axis_exprs, i)
        end
    end

    # Remove duplicate indices on the left hand side
    for i in reverse(eachindex(lhs_indices))
        duplicated = false
        di = lhs_axis_exprs[i]

        # Don't need to check rhs, already done above

        for j = 1:(i - 1)
            if lhs_indices[j] == lhs_indices[i]
                duplicated = true
                dj = lhs_axis_exprs[j]
                push!(dimension_checks, :(@assert $dj == $di))
            end
        end

        if duplicated
            deleteat!(lhs_indices, i)
            deleteat!(lhs_axis_exprs, i)
        end
    end

    # Create output array if specified by user
    @gensym T

    if Meta.isexpr(expr, :(:=))
        rhs_type = :(promote_type($([:(eltype($arr)) for arr in rhs_arrays]...)))
        type_definition = :(local $T = $rhs_type)
        output_definition = if length(lhs_axis_exprs) > 0
            :($(lhs_arrays[1]) = Array{$T}(undef, $(lhs_axis_exprs...)))
        else
            :($(lhs_arrays[1]) = zero($T))
        end
        assignment_op = :(=)
    else
        type_definition = :(local $T = eltype($(lhs_arrays[1])))
        output_definition = :(nothing)
        assignment_op = expr.head
    end

    # Copy the index expression to modify it.
    # loop_expr is the Expr we'll build loops around.
    loop_expr = unquote_offsets!(copy(expr))

    # Nest loops to iterate over the destination variables
    if length(rhs_indices) > 0
        # Innermost expression has form s += rhs
        @gensym s
        loop_expr.args[1] = s
        loop_expr.head = :(+=)

        # Nest loops to iterate over the summed out variables
        loop_expr = nest_loops(loop_expr, rhs_indices, rhs_axis_exprs)

        # Prepend with s = 0, and append with assignment
        # to the left hand side of the equation.
        lhs_assignment = Expr(assignment_op, lhs, s)

        loop_expr = quote
            local $s = zero($T)
            $loop_expr
            $lhs_assignment
        end

        # Now loop over indices appearing on lhs, if any.
        loop_expr = nest_loops(loop_expr, lhs_indices, lhs_axis_exprs)
    else
        # We do not sum over any indices, only loop over lhs.
        loop_expr.head = assignment_op
        loop_expr = nest_loops(loop_expr, lhs_indices, lhs_axis_exprs)
    end

    # Check bounds
    if inbounds
        loop_expr = :(@inbounds $loop_expr)
    end

    full_expression = quote
        $type_definition
        $output_definition
        $(dimension_checks...)

        let $([lhs_indices; rhs_indices]...)
            $loop_expr
        end

        $(lhs_arrays[1])
    end

    return esc(full_expression)
end

# Check validity
function check_index_occurrence(lhs_indices, rhs_indices)
    if !issubset(lhs_indices, rhs_indices)
        missing_indices = setdiff(lhs_indices, rhs_indices)

        if length(missing_indices) == 1
            missing_string = "\"$(missing_indices[1])\""
            throw(ArgumentError(string(
                "Index ", missing_string, " is occuring only on left side")))
        else
            missing_string = join(("\"$ix\"" for ix in missing_indices),
                                  ", ", " and ")
            throw(ArgumentError(string(
                "Indices ", missing_string, " are occuring only on left side")))
        end
    end
end

# Expand loops
function nest_loops(expr::Expr, index_names::Vector{Symbol}, axis_expressions::Vector{Expr})
    isempty(index_names) && return expr

    expr = nest_loop(expr, index_names[1], axis_expressions[1])
    for j = 2:length(index_names)
        expr = nest_loop(expr, index_names[j], axis_expressions[j])
    end

    return expr
end

# Expand loops
function nest_loop(expr::Expr, index_name::Symbol, axis_expression::Expr)
    loop = :(for $index_name = 1:$axis_expression
                 $expr
             end)
    return loop
end

extractindices(expr) = extractindices!(expr, Symbol[], Symbol[], Expr[])

function extractindices!(expr::Symbol,
                         array_names::Vector{Symbol},
                         index_names::Vector{Symbol},
                         axis_expressions::Vector{Expr})
    push!(array_names, expr)
    return array_names, index_names, axis_expressions
end

function extractindices!(expr::Number,
                         array_names::Vector{Symbol},
                         index_names::Vector{Symbol},
                         axis_expressions::Vector{Expr})
    return array_names, index_names, axis_expressions
end

function extractindices!(expr::Expr,
                         array_names::Vector{Symbol},
                         index_names::Vector{Symbol},
                         axis_expressions::Vector{Expr})
    if Meta.isexpr(expr, :ref)
        array_name = expr.args[1]
        push!(array_names, array_name)
        for (dimension, index_expr) in enumerate(expr.args[2:end])
            pushindex!(index_expr,
                       array_name,
                       dimension,
                       array_names,
                       index_names,
                       axis_expressions)
        end
    elseif Meta.isexpr(expr, :call)
        for arg in expr.args[2:end]
            extractindices!(arg, array_names, index_names, axis_expressions)
        end
    else
        throw(ArgumentError("Invalid expression head: `:$(expr.head)`"))
    end

    return return array_names, index_names, axis_expressions
end

function pushindex!(expr::Symbol,
                    array_name::Symbol,
                    dimension::Int,
                    array_names,
                    index_names,
                    axis_expressions)
    push!(index_names, expr)
    push!(axis_expressions, :(size($array_name, $dimension)))
end

function pushindex!(expr::Number,
                    array_name::Symbol,
                    dimension::Int,
                    array_names,
                    index_names,
                    axis_expressions)
    return nothing
end

function pushindex!(expr::Expr,
                    array_name::Symbol,
                    dimension::Int,
                    array_names,
                    index_names,
                    axis_expressions)
    if Meta.isexpr(expr, :call) && length(expr.args) == 3
        op = expr.args[1]

        index_name = expr.args[2]
        @assert typeof(index_name) == Symbol

        offset_expr = expr.args[3]

        if offset_expr isa Integer
            offset = expr.args[3]::Integer
        elseif offset_expr isa Expr && Meta.isexpr(offset_expr, :quote)
            offset = offset_expr.args[1]
        elseif offset_expr isa QuoteNode
            offset = offset_expr.value
        else
            throw(ArgumentError("Improper expression inside reference on rhs"))
        end

        push!(index_names, index_name)

        if op == :+
            push!(axis_expressions, :(size($array_name, $dimension) - $offset))
        elseif op == :-
            push!(axis_expressions, :(size($array_name, $dimension) + $offset))
        else
            throw(ArgumentError(string("Operations inside ref on rhs are ",
                                       "limited to `+` and `-`")))
        end
    elseif Meta.isexpr(expr, :quote)
        return nothing
    else
        throw(ArgumentError("Invalid index expression: `$(expr)`"))
    end
end

function unquote_offsets!(expr::Expr, inside_ref = false)
    inside_ref |= Meta.isexpr(expr, :ref)

    for i in eachindex(expr.args)
        if expr.args[i] isa Expr
            if Meta.isexpr(expr.args[i], :quote) && inside_ref
                expr.args[i] = expr.args[i].args[1]
            else
                unquote_offsets!(expr.args[i], inside_ref)
            end
        end
    end

    return expr
end

#=
### *Math* : *Curve Fitting*
=#

#=
*Remarks* :

The following codes implements the Levenberg-Marquardt algorithm for
curve fitting.

Actually, these codes are borrowed from the following repositories:

* https://github.com/JuliaNLSolvers/LsqFit.jl
* https://github.com/JuliaNLSolvers/OptimBase.jl
* https://github.com/JuliaNLSolvers/NLSolversBase.jl
* https://github.com/JuliaDiff/FiniteDiff.jl

Of cource, these codes are greatly simplified, and only the vital
features are retained. Only the `curve_fit()` function is exported. If
any errors occur, please turn to the original version of `curve_fit()`
as implemented in the `LsqFit.jl` package.
=#

"""
    OnceDifferentiable

Mutable struct. It is used for objectives and solvers where the gradient
is available/exists.

### Members
* ℱ! -> Objective. It is actually a function call and return objective.
* 𝒥! -> It is a function call as well and returns jacobian of objective.
* 𝐹  -> Cache for ℱ! output.
* 𝐽  -> Cache for 𝒥! output.
"""
mutable struct OnceDifferentiable
    ℱ!
    𝒥!
    𝐹
    𝐽
end

"""
    OnceDifferentiable(𝑓, p0::AbstractArray, 𝐹::AbstractArray)

Constructor for OnceDifferentiable struct. `𝑓` is the function, `p0` is
the inital guess, `𝐹 = 𝑓(p0)` is the cache for `𝑓`'s output.
"""
function OnceDifferentiable(𝑓, p0::AbstractArray, 𝐹::AbstractArray)
    # Backup 𝑓(x) to 𝐹.
    function ℱ!(𝐹, x)
        copyto!(𝐹, 𝑓(x))
    end

    # Calculate jacobian for 𝑓(x), the results are stored in 𝐽.
    # The finite difference method is used.
    function 𝒥!(𝐽, x)
        rel_step = cbrt(eps(real(eltype(x))))
        abs_step = rel_step
        @inbounds for i ∈ 1:length(x)
            xₛ = x[i]
            ϵ = max(rel_step * abs(xₛ), abs_step)
            x[i] = xₛ + ϵ
            f₂ = vec(𝑓(x))
            x[i] = xₛ - ϵ
            f₁ = vec(𝑓(x))
            𝐽[:,i] = (f₂ - f₁) ./ (2 * ϵ)
            x[i] = xₛ
        end
    end

    # Create memory space for jacobian matrix
    𝐽 = eltype(p0)(NaN) .* vec(𝐹) .* vec(p0)'

    # Call the default constructor
    OnceDifferentiable(ℱ!, 𝒥!, 𝐹, 𝐽)
end

"""
    value(obj::OnceDifferentiable)

Return `obj.𝐹`. `obj` will not be affected.
"""
value(obj::OnceDifferentiable) = obj.𝐹

"""
    value(obj::OnceDifferentiable, 𝐹, x)

Return `𝑓(x)`. `obj` will not be affected, but `𝐹` is updated.
"""
value(obj::OnceDifferentiable, 𝐹, x) = obj.ℱ!(𝐹, x)

"""
    value!(obj::OnceDifferentiable, x)

Return `𝑓(x)`. `obj.𝐹` will be updated and returned.
"""
function value!(obj::OnceDifferentiable, x)
    obj.ℱ!(obj.𝐹, x)
    obj.𝐹
end

"""
    jacobian(obj::OnceDifferentiable)

Return `obj.𝐽`. `obj` will not be affected.
"""
jacobian(obj::OnceDifferentiable) = obj.𝐽

"""
    jacobian(obj::OnceDifferentiable, 𝐽, x)

Return jacobian. `obj` will not be affected, but `𝐽` is updated.
"""

jacobian(obj::OnceDifferentiable, 𝐽, x) = obj.𝒥!(𝐽, x)

"""
    jacobian!(obj::OnceDifferentiable, x)

Return jacobian. `obj.𝐽` will be updated and returned.
"""
function jacobian!(obj::OnceDifferentiable, x)
    obj.𝒥!(obj.𝐽, x)
    obj.𝐽
end

"""
    LMOptimizationResults{T,N}

It is used to save the optimization results of the `levenberg_marquardt`
algorithm.

### Members
* x₀         -> Initial guess for the solution.
* minimizer  -> Final results for the solution.
* minimum    -> Residual.
* iterations -> Number of iterations.
* xconv      -> If the convergence criterion 1 is satisfied.
* gconv      -> If the convergence criterion 2 is satisfied.
"""
struct LMOptimizationResults{T,N}
    x₀ :: Array{T,N}
    minimizer :: Array{T,N}
    minimum :: T
    iterations :: Int
    xconv :: Bool
    gconv :: Bool
end

"""
    levenberg_marquardt(df::OnceDifferentiable, x₀::AbstractVector{T})

Returns the argmin over x of `sum(f(x).^2)` using the Levenberg-Marquardt
algorithm. The function `f` is encoded in `df`. `x₀` is an initial guess
for the solution.

See also: [`OnceDifferentiable`](@ref).
"""
function levenberg_marquardt(df::OnceDifferentiable, x₀::AbstractVector{T} where T)
    # Some predefined constants
    min_diagonal = 1e-6 # lower bound on values of diagonal matrix
    #
    x_tol   = 1e-08 # search tolerance in x
    g_tol   = 1e-12 # search tolerance in gradient
    maxIter = 1000  # maximum number of iterations
    #
    Λₘ = 1e+16 # minimum trust region radius
    λₘ = 1e-16 # maximum trust region radius
    λ  = eltype(x₀)(10) # (inverse of) initial trust region radius
    λᵢ = 10.0  # λ is multiplied by this factor after step below min quality
    λᵣ = 0.10  # λ is multiplied by this factor after good quality steps
    #
    min_step_quality  = 1e-3 # for steps below this quality, the trust region is shrinked
    good_step_quality = 0.75 # for steps above this quality, the trust region is expanded

    # First evaluation
    # Both df.𝐹 and df.𝐽 are updated.
    # And 𝐹 and 𝐽 become aliases of df.𝐹 and df.𝐽, respectively.
    value!(df, x₀)
    jacobian!(df, x₀)
    𝐹 = value(df)
    𝐽 = jacobian(df)

    # Setup convergence criteria
    converged = false
    xconv = false
    gconv = false
    iter = 0

    # Calculate 𝑓(x₀) and initial residual
    x = copy(x₀)
    trial_f = similar(𝐹)
    C_resid = sum(abs2, 𝐹)

    # Create buffers
    𝐽ᵀ𝐽 = diagm(x)
    𝐽δx = similar(𝐹)

    # Main iteration
    while (~converged && iter < maxIter)
        # Update jacobian 𝐽 for new x
        jacobian!(df, x)

        # Solve the equation: [𝐽ᵀ𝐽 + λ diag(𝐽ᵀ𝐽)] δ = 𝐽ᵀ𝐹
        # What we want to get is δ.
        mul!(𝐽ᵀ𝐽, 𝐽', 𝐽)
        #
        𝐷ᵀ𝐷 = diag(𝐽ᵀ𝐽)
        replace!(x -> x ≤ min_diagonal ? min_diagonal : x, 𝐷ᵀ𝐷)
        #
        @simd for i in eachindex(𝐷ᵀ𝐷)
            @inbounds 𝐽ᵀ𝐽[i,i] += λ * 𝐷ᵀ𝐷[i]
        end
        #
        δx = - 𝐽ᵀ𝐽 \ (𝐽' * 𝐹)

        # If the linear assumption is valid, the new residual is predicted.
        mul!(𝐽δx, 𝐽, δx)
        𝐽δx .= 𝐽δx .+ 𝐹
        P_resid = sum(abs2, 𝐽δx)

        # Try to calculate new x, and then 𝐹 ≡ 𝑓(x), and then the residual.
        xnew = x + δx
        value(df, trial_f, xnew)
        T_resid = sum(abs2, trial_f)

        # Step quality = residual change / predicted residual change
        ρ = (T_resid - C_resid) / (P_resid - C_resid)
        if ρ > min_step_quality
            # Update x, 𝑓(x), and residual.
            x .= xnew
            value!(df, x)
            C_resid = T_resid

            # Increase trust region radius
            if ρ > good_step_quality
                λ = max(λᵣ * λ, λₘ)
            end
        else
            # Decrease trust region radius
            λ = min(λᵢ * λ, Λₘ)
        end

        # Increase the iteration
        iter += 1

        # Check convergence criteria:
        # 1. Small gradient: norm(𝐽ᵀ * 𝐹, Inf) < g_tol
        if norm(𝐽' * 𝐹, Inf) < g_tol
            gconv = true
        end
        # 2. Small step size: norm(δx) < x_tol
        if norm(δx) < x_tol * (x_tol + norm(x))
            xconv = true
        end
        # 3. Calculate converged
        converged = gconv | xconv
    end

    # Return the results
    LMOptimizationResults(
        x₀,      # x₀
        x,       # minimizer
        C_resid, # residual
        iter,    # iterations
        xconv,   # xconv
        gconv,   # gconv
    )
end

"""
    LsqFitResult

It encapsulates the results for curve fitting.

### Members
* param     -> Fitted results, i.e, the fitting parameters.
* resid     -> Residuals.
* jacobian  -> Jacobian matrix.
* converged -> If the curve-fitting algorithm is converged.
"""
struct LsqFitResult
    param
    resid
    jacobian
    converged
end

"""
    curve_fit(model, x, y, p0)

Fit data to a non-linear `model`. `p0` is an initial model parameter guess.
The return object is a composite type (`LsqFitResult`).

See also: [`LsqFitResult`](@ref).
"""
function curve_fit(model, x::AbstractArray, y::AbstractArray, p0::AbstractArray)
    f = (p) -> model(x, p) - y
    r = f(p0)
    R = OnceDifferentiable(f, p0, r)
    OR = levenberg_marquardt(R, p0)
    p = OR.minimizer
    conv = OR.xconv || OR.gconv
    return LsqFitResult(p, value!(R, p), jacobian!(R, p), conv)
end

export OnceDifferentiable
export LMOptimizationResults
export LsqFitResult
export value
export value!
export jacobian
export jacobian!
export levenberg_marquardt

#=
### *Math* : *Numerical Optimization*
=#

#=
*Remarks* :

The following codes implements the BFGS algorithm for numerical optimization.

Actually, these codes are borrowed from the following repositories:

* https://github.com/JuliaNLSolvers/Optim.jl
* https://github.com/JuliaNLSolvers/NLSolversBase.jl
* https://github.com/JuliaNLSolvers/LineSearches.jl

Of cource, these codes are greatly simplified, and only the vital
features are retained. Only the `optimize()` function is exported. If
any errors occur, please turn to the original version of `optimize()`
as implemented in the `Optim.jl` package.
=#

"""
    BFGSDifferentiable

Mutable struct. It is used for objectives and solvers where the gradient
is available/exists.

### Members
* ℱ! -> Objective. It is actually a function call and return objective.
* 𝒟! -> It is a function call as well and returns derivative of objective.
* 𝐹  -> Cache for ℱ! output.
* 𝐷  -> Cache for 𝒟! output.
"""
mutable struct BFGSDifferentiable
    ℱ!
    𝒟!
    𝐹
    𝐷
end

"""
    BFGSDifferentiable(f, df, x::AbstractArray)

Constructor for BFGSDifferentiable struct. `f` is the function, `df` is
the derivative of objective, `x` is the initial guess.
"""
function BFGSDifferentiable(f, df, x::AbstractArray)
    𝐹 = real(zero(eltype(x)))
    T = promote_type(eltype(x), eltype(𝐹))
    𝐷 = fill!(T.(x), T(NaN))
    return BFGSDifferentiable(f, df, copy(𝐹), copy(𝐷))
end

"""
    value(obj::BFGSDifferentiable)

Return `obj.𝐹`. `obj` will not be affected.
"""
value(obj::BFGSDifferentiable) = obj.𝐹

"""
    gradient(obj::BFGSDifferentiable)

Return `obj.𝐷`. `obj` will not be affected.
"""
gradient(obj::BFGSDifferentiable) = obj.𝐷

"""
    value_gradient!(obj::BFGSDifferentiable, x)

Evaluate objective and derivative at `x`. `obj.𝐹` and `obj.𝐷` should be
updated. Note that here `obj.𝒟!` is actually `nac.jl/smooth_norm()`.
"""
function value_gradient!(obj::BFGSDifferentiable, x)
    # Note that gradient(obj), i.e obj.𝐷, should be updated in obj.𝒟!().
    obj.𝒟!(gradient(obj), x)
    obj.𝐹 = obj.ℱ!(x)
end

"""
    BFGSState

Mutable struct. It is used to trace the history of states visited.

### Members
* x     -> Current position.
* ls    -> Current search direction.
* δx    -> Changes in position.
* δg    -> Changes in gradient.
* xₚ    -> Previous position.
* gₚ    -> Previous gradient.
* fₚ    -> Previous f (f in xₚ).
* H⁻¹   -> Current inverse Hessian matrix.
* alpha -> A internal parameter to control the BFGS algorithm.
"""
mutable struct BFGSState{Tx, Tm, T, G}
    x :: Tx
    ls :: Tx
    δx :: Tx
    δg :: Tx
    xₚ :: Tx
    gₚ :: G
    fₚ :: T
    H⁻¹ :: Tm
    alpha :: T
end

"""
    BFGSOptimizationResults

It is used to save the optimization results of the BFGS algorithm.

### Members
* x₀         -> Initial guess for the solution.
* minimizer  -> Final results for the solution.
* minimum    -> Objective at the final solution.
* iterations -> Number of iterations.
* δx         -> Absolute change in x.
* Δx         -> Relative change in x.
* δf         -> Absolute change in f.
* Δf         -> Relative change in f.
* resid      -> Maximum gradient of f at the final solution.
* gconv      -> If the convergence criterion is satisfied
"""
mutable struct BFGSOptimizationResults{Tx, Tc, Tf}
    x₀ :: Tx
    minimizer  :: Tx
    minimum    :: Tf
    iterations :: Int
    δx :: Tc
    Δx :: Tc
    δf :: Tc
    Δf :: Tc
    resid :: Tc
    gconv :: Bool
end

"""
    maxdiff(x::AbstractArray, y::AbstractArray)

Return the maximum difference between two arrays: `x` and `y`. Note that
the sizes of `x` and `y` should be the same.
"""
function maxdiff(x::AbstractArray, y::AbstractArray)
    return mapreduce((a, b) -> abs(a - b), max, x, y)
end

"""
    eval_δf(d::BFGSDifferentiable, s::BFGSState)

Evaluate the absolute changes in f.
"""
eval_δf(d::BFGSDifferentiable, s::BFGSState) = abs(value(d) - s.fₚ)

"""
    eval_Δf(d::BFGSDifferentiable, s::BFGSState)

Evaluate the relative changes in f.
"""
eval_Δf(d::BFGSDifferentiable, s::BFGSState) = eval_δf(d, s) / abs(value(d))

"""
    eval_δx(s::BFGSState)

Evaluate the absolute changes in x.
"""
eval_δx(s::BFGSState) = maxdiff(s.x, s.xₚ)

"""
    eval_Δx(s::BFGSState)

Evaluate the relative changes in x.
"""
eval_Δx(s::BFGSState) = eval_δx(s) / maximum(abs, s.x)

"""
    eval_resid(d::BFGSDifferentiable)

Evaluate residual (maximum gradient of f at the current position).
"""
eval_resid(d::BFGSDifferentiable) = maximum(abs, gradient(d))

"""
    optimize(f, g, x₀::AbstractArray; max_iter::I64 = 1000)

Return the argmin over x of `f(x)` using the BFGS algorithm. Here, `f`
is the function call, and `g` will return the gradient of `f`, `x₀` is
an initial guess for the solution.
"""
function optimize(f, g, x₀::AbstractArray; max_iter::I64 = 1000)
    # Initialize time stamp
    t₀ = time()

    # Create BFGSDifferentiable
    d = BFGSDifferentiable(f, g, x₀)

    # Create BFGSState
    s = init_state(d, x₀)

    # Prepare counter
    iteration = 0

    # Print trace for optimization progress
    @printf("Tracing BFGS Optimization\n")
    trace(d, iteration, time() - t₀)

    # Setup convergence flag
    gconv = !isfinite(value(d)) || any(!isfinite, gradient(d))
    while !gconv && iteration < max_iter
        iteration += 1

        # Update line search direction
        ls_success = !update_state!(d, s)
        if !ls_success
            break
        end

        # Update the value of f and its gradient
        update_g!(d, s)

        # Update the Hessian matrix
        update_h!(d, s)

        # Print trace for optimization progress
        trace(d, iteration, time() - t₀)

        # Check the gradient
        if !all(isfinite, gradient(d))
            @warn "Terminated early due to NaN in gradient."
            break
        end

        # Check whether convergence criterion is satisfied
        gconv = (eval_resid(d) ≤ 1e-8)
    end

    # Return BFGSOptimizationResults
    BFGSOptimizationResults(x₀, s.x, value(d), iteration,
                            eval_δx(s), eval_Δx(s),
                            eval_δf(d, s), eval_Δf(d, s),
                            eval_resid(d),
                            gconv)
end

"""
    init_state(d::BFGSDifferentiable, x₀::AbstractArray)

Create a BFGSState object. Note that `d` should be updated in this
function (`d.𝐹` and `d.𝐷`). `x₀` is an initial guess for the solution.

See also: [`BFGSDifferentiable`](@ref), [`BFGSState`](@ref).
"""
function init_state(d::BFGSDifferentiable, x₀::AbstractArray)
    # Update `d.𝐹` and `d.𝐷` using x₀.
    value_gradient!(d, x₀)

    # Prepare inverse Hessian matrix
    T = eltype(x₀)
    x_ = reshape(x₀, :)
    H⁻¹ = x_ .* x_' .* false
    idxs = diagind(H⁻¹)
    scale = T(1)
    @. @view(H⁻¹[idxs]) = scale * true

    # Return BFGSState
    BFGSState(x₀, similar(x₀), similar(x₀), similar(x₀), copy(x₀),
              copy(gradient(d)), real(T)(NaN),
              H⁻¹,
              real(one(T)))
end

"""
    update_state!(d::BFGSDifferentiable, s::BFGSState)

Evaluate line search direction and change of x. New position and old
gradient are saved in `s`.

See also: [`BFGSDifferentiable`](@ref), [`BFGSState`](@ref).
"""
function update_state!(d::BFGSDifferentiable, s::BFGSState)
    T = eltype(s.ls)

    # Set the search direction
    #
    # Note that search direction is the negative gradient divided by
    # the approximate Hessian matrix.
    mul!(vec(s.ls), s.H⁻¹, vec(gradient(d)))
    rmul!(s.ls, T(-1))

    # Maintain a record of the previous gradient
    copyto!(s.gₚ, gradient(d))

    # Determine the distance of movement along the search line
    lssuccess = linesearch!(d, s)

    # Update current position
    s.δx .= s.alpha .* s.ls
    s.x .= s.x .+ s.δx

    # Break on linesearch error
    lssuccess == false
end

"""
    update_g!(d::BFGSDifferentiable, s::BFGSState)

Update the function value and gradient (`d.𝐹` and `d.𝐷` are changed).

See also: [`BFGSDifferentiable`](@ref), [`BFGSState`](@ref).
"""
update_g!(d::BFGSDifferentiable, s::BFGSState) = value_gradient!(d, s.x)

"""
    update_h!(d::BFGSDifferentiable, s::BFGSState)

Try to evaluate the new Hessian matrix.

See also: [`BFGSDifferentiable`](@ref), [`BFGSState`](@ref).
"""
function update_h!(d::BFGSDifferentiable, s::BFGSState)
    n = length(s.x)
    su = similar(s.x)

    # Measure the change in the gradient
    s.δg .= gradient(d) .- s.gₚ

    # Update the inverse Hessian approximation by using the
    # famous Sherman-Morrison equation
    dx_dg = real(dot(s.δx, s.δg))
    if dx_dg > 0
        mul!(vec(su), s.H⁻¹, vec(s.δg))

        c1 = (dx_dg + real(dot(s.δg, su))) / (dx_dg' * dx_dg)
        c2 = 1 / dx_dg

        # H⁻¹ = H⁻¹ + c1 * (s * s') - c2 * (u * s' + s * u')
        if s.H⁻¹ isa Array
            H⁻¹ = s.H⁻¹; dx = s.δx; u = su;
            @inbounds for j in 1:n
                c1dxj = c1 * dx[j]'
                c2dxj = c2 * dx[j]'
                c2uj  = c2 *  u[j]'
                for i in 1:n
                    H⁻¹[i, j] = muladd(dx[i], c1dxj,
                                    muladd(-u[i], c2dxj,
                                        muladd(c2uj, -dx[i], H⁻¹[i, j])))
                end
            end
        else
            mul!(s.H⁻¹, vec(s.δx), vec(s.δx)', +c1, 1)
            mul!(s.H⁻¹, vec(su  ), vec(s.δx)', -c2, 1)
            mul!(s.H⁻¹, vec(s.δx), vec(su  )', -c2, 1)
        end
    end
end

"""
    trace(d::BFGSDifferentiable, iter::I64, curr_time::F64)

Print some useful information about the optimization procedure.
"""
function trace(d::BFGSDifferentiable, iter::I64, curr_time::F64)
    gnorm = norm(gradient(d), Inf)
    #
    @printf("iter = %4d ", iter)
    @printf("𝑓 = %8.4e ", value(d))
    @printf("||∂𝑓/∂x|| = %8.4e ", gnorm)
    @printf("time = %8.4f (s)\n", curr_time)
    #
    flush(stdout)
end

"""
    linesearch!(d::BFGSDifferentiable, s::BFGSState)

Evaluate line search direction. Actually, `s.alpha`, `s.fₚ`, and `s.xₚ`
will be updated in this function.

See also: [`BFGSDifferentiable`](@ref), [`BFGSState`](@ref).
"""
function linesearch!(d::BFGSDifferentiable, s::BFGSState)
    # Calculate search direction dϕ₀
    dϕ₀ = real(dot(gradient(d), s.ls))

    # Reset the direction if it becomes corrupted
    if dϕ₀ >= zero(dϕ₀)
        dϕ₀ = real(dot(gradient(d), s.ls))
    end

    # Guess an alpha
    LS(s, 1.0, false)

    # Store current x and f(x) for next iteration
    ϕ₀ = value(d)
    s.fₚ = ϕ₀
    copyto!(s.xₚ, s.x)

    # Perform line search using the Hager-Zhang algorithm
    try
        s.alpha, _ = LS(d, s.x, s.ls, s.alpha, ϕ₀, dϕ₀)
        return true # lssuccess = true
    catch ex
        # Catch LineSearchException to allow graceful exit
        if isa(ex, LineSearchException)
            s.alpha = ex.alpha
            return false # lssuccess = false
        else
            rethrow(ex)
        end
    end
end

"""
    converged(r::BFGSOptimizationResults)

Check whether the optimization is converged.

See also: [`BFGSOptimizationResults`](@ref).
"""
function converged(r::BFGSOptimizationResults)
    conv_flags = r.gconv
    x_isfinite = isfinite(r.δx) || isnan(r.Δx)
    f_isfinite = if r.iterations > 0
            isfinite(r.δf) || isnan(r.Δf)
        else
            true
        end
    g_isfinite = isfinite(r.resid)
    return conv_flags && all((x_isfinite, f_isfinite, g_isfinite))
end

export BFGSDifferentiable
export BFGSState
export BFGSOptimizationResults
export value
export gradient
export value_gradient!
export maxdiff
export eval_δf
export eval_Δf
export eval_δx
export eval_Δx
export eval_resid
export init_state
export update_state!
export update_g!
export update_h!
export trace
export linesearch!
export converged

#=
### *Math* : *Line Search*
=#

#
# Conjugate gradient line search implementation from:
#   W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a
#     conjugate gradient method with guaranteed descent. ACM
#     Transactions on Mathematical Software 32: 113–137.
#

#
# The following codes are borrowed from:
#     https://github.com/JuliaNLSolvers/LineSearches.jl
#

"""
    LineSearchException

Mutable struct. It contains information about the error occured in
the line search.

### Members
* message -> Error message.
* alpha   -> A key parameter used to control line search.
"""
mutable struct LineSearchException{T<:Real} <: Exception
    message :: AbstractString
    alpha :: T
end

"""
    LS(state::BFGSState, alpha::T, scaled::Bool)

Line search: initial and static version.
"""
function LS(state::BFGSState, alpha::T, scaled::Bool) where T
    PT = promote_type(T, real(eltype(state.ls)))
    if scaled == true && (ns = real(norm(state.ls))) > convert(PT, 0)
        state.alpha = convert(PT, min(alpha, ns)) / ns
    else
        state.alpha = convert(PT, alpha)
    end
end

"""
    LS(df::BFGSDifferentiable,
        x::Vector{C64}, s::Vector{C64},
        c::F64, phi_0::F64, dphi_0::F64)

Line search: Hager-Zhang algorithm.
"""
function LS(df::BFGSDifferentiable,
            x::Vector{C64}, s::Vector{C64},
            c::F64, phi_0::F64, dphi_0::F64)
    delta = 0.1
    sigma = 0.9
    alphamax = Inf
    rho = 5.0
    epsilon = 1e-6
    gamma = 0.66
    linesearchmax = 50
    psi3 = 0.1
    mayterminate = Ref{Bool}(false)

    ϕdϕ = make_ϕdϕ(df, similar(x), x, s)

    T = typeof(c)
    zeroT = convert(T, 0)
    #
    if !(isfinite(phi_0) && isfinite(dphi_0))
        throw(LineSearchException("Value and slope at step length = 0 must be finite.", T(0)))
    end
    #
    if dphi_0 >= eps(T) * abs(phi_0)
        throw(LineSearchException("Search direction is not a direction of descent.", T(0)))
    elseif dphi_0 >= 0
        return zeroT, phi_0
    end

    # Prevent values of x_new = x+αs that are likely to
    # make ϕ(x_new) infinite
    iterfinitemax::Int = ceil(Int, -log2(eps(T)))
    alphas = [zeroT] # for bisection
    values = [phi_0]
    slopes = [dphi_0]

    phi_lim = phi_0 + epsilon * abs(phi_0)
    @assert c >= 0
    c <= eps(T) && return zeroT, phi_0
    @assert isfinite(c) && c <= alphamax
    phi_c, dphi_c = ϕdϕ(c)
    iterfinite = 1
    #
    while !(isfinite(phi_c) && isfinite(dphi_c)) && iterfinite < iterfinitemax
        mayterminate[] = false
        iterfinite += 1
        c *= psi3
        phi_c, dphi_c = ϕdϕ(c)
    end
    #
    if !(isfinite(phi_c) && isfinite(dphi_c))
        @warn("Failed to achieve finite new evaluation point, using alpha=0")
        mayterminate[] = false # reset in case another initial guess is used next
        return zeroT, phi_0
    end
    #
    push!(alphas, c)
    push!(values, phi_c)
    push!(slopes, dphi_c)

    # If c was generated by quadratic interpolation, check whether it
    # satisfies the Wolfe conditions
    if mayterminate[] &&
          satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma)
        # Reset in case another initial guess is used next
        mayterminate[] = false
        return c, phi_c # phi_c
    end

    # Initial bracketing step (HZ, stages B0-B3)
    isbracketed = false
    ia = 1
    ib = 2
    @assert length(alphas) == 2
    iter = 1
    cold = -one(T)
    while !isbracketed && iter < linesearchmax
        if dphi_c >= zeroT
            # We've reached the upward slope, so we have b; examine
            # previous values to find a
            ib = length(alphas)
            for i = (ib - 1):-1:1
                if values[i] <= phi_lim
                    ia = i
                    break
                end
            end
            isbracketed = true
        elseif values[end] > phi_lim
            # The value is higher, but the slope is downward, so we must
            # have crested over the peak. Use bisection.
            ib = length(alphas)
            ia = 1
            if c ≉  alphas[ib] || slopes[ib] >= zeroT
                error("c = ", c)
            end
            ia, ib = ls_bisect!(ϕdϕ, alphas, values, slopes, ia, ib, phi_lim)
            isbracketed = true
        else
            # We'll still going downhill, expand the interval and try again.
            # Reaching this branch means that dphi_c < 0 and phi_c <= phi_0 + ϵ_k
            # So cold = c has a lower objective than phi_0 up to epsilon.
            # This makes it a viable step to return if bracketing fails.

            # Bracketing can fail if no cold < c <= alphamax can be found
            # with finite phi_c and dphi_c. Going back to the loop with
            # c = cold will only result in infinite cycling. So returning
            # (cold, phi_cold) and exiting the line search is the best move.
            cold = c
            phi_cold = phi_c
            if nextfloat(cold) >= alphamax
                # Reset in case another initial guess is used next
                mayterminate[] = false
                return cold, phi_cold
            end
            c *= rho
            if c > alphamax
                c = alphamax
            end
            phi_c, dphi_c = ϕdϕ(c)
            iterfinite = 1
            while !(isfinite(phi_c) && isfinite(dphi_c)) &&
                    c > nextfloat(cold) && iterfinite < iterfinitemax
                # Shrinks alphamax, assumes that steps >= c can never
                # have finite phi_c and dphi_c.
                alphamax = c
                iterfinite += 1
                c = (cold + c) / 2
                phi_c, dphi_c = ϕdϕ(c)
            end
            if !(isfinite(phi_c) && isfinite(dphi_c))
                return cold, phi_cold
            end
            push!(alphas, c)
            push!(values, phi_c)
            push!(slopes, dphi_c)
        end
        iter += 1
    end

    while iter < linesearchmax
        a = alphas[ia]
        b = alphas[ib]
        @assert b > a
        if b - a <= eps(b)
            # Reset in case another initial guess is used next
            mayterminate[] = false
            return a, values[ia] # lsr.value[ia]
        end
        iswolfe, iA, iB = ls_secant2!(ϕdϕ, alphas, values, slopes, ia, ib, phi_lim, delta, sigma)
        if iswolfe
            # Reset in case another initial guess is used next
            mayterminate[] = false
            return alphas[iA], values[iA]
        end
        A = alphas[iA]
        B = alphas[iB]
        @assert B > A
        if B - A < gamma * (b - a)
            if nextfloat(values[ia]) >= values[ib] && nextfloat(values[iA]) >= values[iB]
                # It's so flat, secant didn't do anything useful, time to quit
                # Reset in case another initial guess is used next
                mayterminate[] = false
                return A, values[iA]
            end
            ia = iA
            ib = iB
        else
            # Secant is converging too slowly, use bisection
            c = (A + B) / convert(T, 2)

            phi_c, dphi_c = ϕdϕ(c)
            @assert isfinite(phi_c) && isfinite(dphi_c)
            push!(alphas, c)
            push!(values, phi_c)
            push!(slopes, dphi_c)

            ia, ib = ls_update!(ϕdϕ, alphas, values, slopes, iA, iB, length(alphas), phi_lim)
        end
        iter += 1
    end

    throw(LineSearchException("Linesearch failed to converge, reached maximum
                               iterations $(linesearchmax).", alphas[ia]))
end

function make_ϕdϕ(df::BFGSDifferentiable, x_new, x, s)
    function ϕdϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate ∇f(x+α*s)
        value_gradient!(df, x_new)

        # Calculate ϕ'(a_i)
        value(df), real(dot(gradient(df), s))
    end
    ϕdϕ
end

# Check Wolfe & approximate Wolfe
function satisfies_wolfe(c::F64, phi_c::F64, dphi_c::F64,
                         phi_0::F64, dphi_0::F64, phi_lim::F64,
                         delta::F64, sigma::F64)
    wolfe1 = delta * dphi_0 >= (phi_c - phi_0) / c &&
               dphi_c >= sigma * dphi_0
    wolfe2 = (2 * delta - 1) * dphi_0 >= dphi_c >= sigma * dphi_0 &&
               phi_c <= phi_lim
    return wolfe1 || wolfe2
end

# HZ, stage U3 (with theta=0.5)
function ls_bisect!(ϕdϕ, alphas::Vector{F64}, values::Vector{F64},
                    slopes::Vector{F64}, ia::I64, ib::I64, phi_lim::F64)
    T = eltype(alphas)
    gphi = convert(T, NaN)
    a = alphas[ia]
    b = alphas[ib]

    # Debugging (HZ, conditions shown following U3)
    zeroT = convert(T, 0)
    #
    @assert slopes[ia] < zeroT
    @assert values[ia] <= phi_lim
    @assert slopes[ib] < zeroT
    @assert values[ib] > phi_lim
    @assert b > a
    #
    while b - a > eps(b)
        d = (a + b) / convert(T, 2)
        phi_d, gphi = ϕdϕ(d)
        @assert isfinite(phi_d) && isfinite(gphi)

        push!(alphas, d)
        push!(values, phi_d)
        push!(slopes, gphi)

        id = length(alphas)
        #
        if gphi >= zeroT
            return ia, id # replace b, return
        end
        #
        if phi_d <= phi_lim
            a = d # replace a, but keep bisecting until dphi_b > 0
            ia = id
        else
            b = d
            ib = id
        end
    end

    return ia, ib
end

# HZ, stages S1-S4
function ls_secant(a::F64, b::F64, dphi_a::F64, dphi_b::F64)
    return (a * dphi_b - b * dphi_a) / (dphi_b - dphi_a)
end

function ls_secant2!(ϕdϕ, alphas::Vector{F64},
                     values::Vector{F64}, slopes::Vector{F64},
                     ia::I64, ib::I64,
                     phi_lim::F64, delta::F64, sigma::F64)
    phi_0 = values[1]
    dphi_0 = slopes[1]
    a = alphas[ia]
    b = alphas[ib]
    dphi_a = slopes[ia]
    dphi_b = slopes[ib]

    T = eltype(slopes)
    zeroT = convert(T, 0)
    #
    if !(dphi_a < zeroT && dphi_b >= zeroT)
        error(string("Search direction is not a direction of descent; ",
                     "this error may indicate that user-provided derivatives are inaccurate. ",
                      @sprintf "(dphi_a = %f; dphi_b = %f)" dphi_a dphi_b))
    end
    #
    c = ls_secant(a, b, dphi_a, dphi_b)
    @assert isfinite(c)
    #
    phi_c, dphi_c = ϕdϕ(c)
    @assert isfinite(phi_c) && isfinite(dphi_c)

    push!(alphas, c)
    push!(values, phi_c)
    push!(slopes, dphi_c)

    ic = length(alphas)
    if satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma)
        return true, ic, ic
    end

    iA, iB = ls_update!(ϕdϕ, alphas, values, slopes, ia, ib, ic, phi_lim)
    a = alphas[iA]
    b = alphas[iB]

    if iB == ic
        # We updated b, make sure we also update a
        c = ls_secant(alphas[ib], alphas[iB], slopes[ib], slopes[iB])
    elseif iA == ic
        # We updated a, do it for b too
        c = ls_secant(alphas[ia], alphas[iA], slopes[ia], slopes[iA])
    end
    #
    if (iA == ic || iB == ic) && a <= c <= b
        phi_c, dphi_c = ϕdϕ(c)
        @assert isfinite(phi_c) && isfinite(dphi_c)

        push!(alphas, c)
        push!(values, phi_c)
        push!(slopes, dphi_c)

        ic = length(alphas)
        # Check arguments here
        if satisfies_wolfe(c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, delta, sigma)
            return true, ic, ic
        end
        iA, iB = ls_update!(ϕdϕ, alphas, values, slopes, iA, iB, ic, phi_lim)
    end

    return false, iA, iB
end

# HZ, stages U0-U3
#
# Given a third point, pick the best two that retain the bracket
# around the minimum (as defined by HZ, eq. 29)
# b will be the upper bound, and a the lower bound
function ls_update!(ϕdϕ, alphas::Vector{F64},
                    values::Vector{F64}, slopes::Vector{F64},
                    ia::I64, ib::I64, ic::I64, phi_lim::F64)
    a = alphas[ia]
    b = alphas[ib]
    T = eltype(slopes)
    zeroT = convert(T, 0)

    # Debugging (HZ, eq. 4.4):
    @assert slopes[ia] < zeroT
    @assert values[ia] <= phi_lim
    @assert slopes[ib] >= zeroT
    @assert b > a
    #
    c = alphas[ic]
    phi_c = values[ic]
    dphi_c = slopes[ic]
    #
    if c < a || c > b
        return ia, ib #, 0, 0  # it's out of the bracketing interval
    end
    #
    if dphi_c >= zeroT
        return ia, ic #, 0, 0  # replace b with a closer point
    end

    # We know dphi_c < 0. However, phi may not be monotonic between a
    # and c, so check that the value is also smaller than phi_0.  (It's
    # more dangerous to replace a than b, since we're leaving the
    # secure environment of alpha=0; that's why we didn't check this
    # above.)
    if phi_c <= phi_lim
        return ic, ib#, 0, 0  # replace a
    end

    # phi_c is bigger than phi_0, which implies that the minimum
    # lies between a and c. Find it via bisection.
    return ls_bisect!(ϕdϕ, alphas, values, slopes, ia, ic, phi_lim)
end

export LineSearchException
export LS
