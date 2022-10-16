#
# Project : Gardenia
# Source  : math.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/10/16
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
    newton(fun::Function, guess, kwargs...;
           maxiter::I64 = 20000, mixing::F64 = 0.5)

It implements the well-known newton algorithm to locate root of a given
polynomial function. Here, `fun` means the function, `guess` is the initial
solution, and `kwargs...` denotes the required arguments for `fun`. Please
be careful, `func` is a multiple variable function. It not only returns
the value, but also the jacobian matrix of the function.

See also: [`secant`](@ref).
"""
function newton(fun::Function, guess, kwargs...;
                maxiter::I64 = 20000, mixing::F64 = 0.5)
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
    trapz(x::AbstractMesh, y::AbstractVector{T})

Perform numerical integration by using the composite trapezoidal rule.

See also: [`simpson`](@ref).
"""
function trapz(x::AbstractMesh, y::AbstractVector{T} where T)
    value = dot(x.weight, y)
    return value
end

"""
    trapz(x::Vector{F64}, y::Vector{T}, linear::Bool = false)

Perform numerical integration by using the composite trapezoidal rule.

See also: [`simpson`](@ref).
"""
function trapz(x::Vector{F64}, y::Vector{T} where T, linear::Bool = false)
    # For linear mesh
    if linear
        h = x[2] - x[1]
        value = y[1] + y[end] + 2.0 * sum(y[2:end-1])
        value = h * value / 2.0
    # For non-equidistant mesh
    else
        len = length(x)
        value = 0.0
        for i = 1:len-1
            value = value + (y[i] + y[i+1]) * (x[i+1] - x[i])
        end
        value = value / 2.0
    end

    return value
end

"""
    simpson(x::AbstractVector{F64}, y::AbstractVector{T})

Perform numerical integration by using the simpson rule. Note that the
length of `x` and `y` must be odd numbers. And `x` must be a linear and
uniform mesh.

See also: [`trapz`](@ref).
"""
function simpson(x::AbstractVector{F64}, y::AbstractVector{T} where T)
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
    OptimizationResults{T,N}

It is used to save the optimization results of the levenberg_marquardt
algorithm.

### Members

* x₀         -> Initial guess for the solution.
* minimizer  -> Final results for the solution.
* minimum    -> Residual.
* iterations -> Number of iterations.
* xconv      -> If the convergence criterion 1 is satisfied.
* gconv      -> If the convergence criterion 2 is satisfied.
"""
struct OptimizationResults{T,N}
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
    OptimizationResults(
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
export OptimizationResults
export LsqFitResult
export levenberg_marquardt
export value
export value!
export jacobian
export jacobian!
