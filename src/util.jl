#
# Project : Gardenia
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/09
#

#=
### *Basic Macros*
=#

"""
    @cswitch(constexpr, body)

Provides a C-like switch statement with the *falling through* behavior.
This implementation was borrowed from the following github repository:

* https://github.com/Gnimuc/CSyntax.jl

### Examples
```julia
engine = get_d("engine")
@cswitch engine begin
    @case "vasp"
        just_do_it()
        break

    @default
        sorry()
        break
end
```
"""
macro cswitch(constexpr, body)
    case2label = Dict{Any,Symbol}()
    flow = Expr(:block)
    end_label = gensym("end")
    default_label = end_label

    for arg in body.args
        if Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@case")
            label = gensym("case")
            case2label[arg.args[3]] = label
            labelexpr = Expr(:symboliclabel, label)
            push!(flow.args, labelexpr)
        elseif Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@default")
            default_label = gensym("default")
            labelexpr = Expr(:symboliclabel, default_label)
            push!(flow.args, labelexpr)
        elseif arg == Expr(:break)
            labelexpr = Expr(:symbolicgoto, end_label)
            push!(flow.args, labelexpr)
        else
            push!(flow.args, arg)
        end
    end
    push!(flow.args, Expr(:symboliclabel, end_label))

    jumptable = Expr(:block)
    for (case, label) in case2label
        condition = Expr(:call, :(==), constexpr, case)
        push!(jumptable.args, Expr(:if, condition, Expr(:symbolicgoto, label)))
    end
    push!(jumptable.args[end].args, Expr(:symbolicgoto, default_label))

    return esc(Expr(:block, jumptable, flow))
end

"""
    @time_call(ex)

Evaluate a function call (`ex`), and then print the elapsed time (number
of seconds) it took to execute.

This macro is a variation of the standard `@elapsed` macro.
"""
macro time_call(ex)
    quote
        while false; end
        local t₀ = time_ns()
        $(esc(ex))
        δt = (time_ns() - t₀) / 1e9
        println("Report: Total elapsed time $(δt) s\n")
        flush(stdout)
    end
end

"""
    @pcs(x...)

Try to print colorful strings. Here `x` is a combination of strings and
colors. Its format likes `string1 color1 string2 color2 (repeat)`. For
the supported colors, please check the global dict `COLORS`.

### Examples
```julia-repl
julia> @pcs "Hello world!" blue
julia> @pcs "Hello " red "world!" green
```

See also: [`COLORS`](@ref), [`welcome`](@ref).
"""
macro pcs(x...)
    ex = quote
        # The `args` is actually a Tuple
        args = $x

        # We have to make sure the strings and colors are paired.
        @assert iseven(length(args))

        for i = 1:2:length(args)
            # Construct and check string
            # Sometimes args[i] contains interpolated variables, its
            # type is `Expr`. At this time, we have to evaluate this
            # `Expr` at first to convert it to a format `String`.
            str   = eval(args[i])
            @assert str isa AbstractString
            #
            # Construct and check color
            color = args[i+1]
            @assert color isa Symbol

            # Generate expression
            print(eval(color)(str))
        end
    end

    return :( $(esc(ex)) )
end

#=
### *Query Runtime Environment*
=#

"""
    require()

Check the version of julia runtime environment. It should be higher
than v1.6.x. One of the most important philosophies of the `ACFlow`
package is minimizing the dependence on the third-party libraries as
far as possible. Note that the `ACFlow` package relys on the `TOML`
package to parse the *.toml file. Only in v1.6.0 and higher versions,
julia includes the `TOML` package in its standard library.
"""
function require()
    if VERSION < v"1.6-"
        error("Please upgrade your julia to v1.6.0 or higher")
    end
end

"""
    setup_args(x::Vararg{String})

Setup `ARGS` manually. This function is used only in `REPL` environment.
We can use this function to update `ARGS`, so that the `query_args()`
and the other related functions can work correctly.

### Examples
```julia-repl
julia> setup_args("ac.toml")
1-element Array{String,1}:
 "SrVO3.toml"
```

See also: [`query_args`](@ref).
"""
function setup_args(x::Vararg{String})
    # Make sure it is the REPL
    @assert isinteractive()

    # Clean `ARGS`
    empty!(ARGS)

    # Convert the arguments to an array of strings
    X = collect(x)

    # Push them into `ARGS` one by one
    for i in eachindex(X)
        push!(ARGS, X[i])
    end

    # Return `ARGS`, only for debug.
    ARGS
end

"""
    query_args()

Check whether the configuration file (`case.toml`) is provided.

See also: [`setup_args`](@ref).
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
        error("Please specify the configuration file")
    else
        ARGS[1]
    end
end

#=
### *Colorful Outputs*
=#

"""
    welcome()

Print out the welcome messages to the screen.
"""
function welcome()
    println(  red("╔═╗╔═╗╔═╗"), magenta("┬  ┌─┐┬ ┬"))
    println(green("╠═╣║  ╠╣ "), magenta("│  │ ││││"))
    println( blue("╩ ╩╚═╝╚  "), magenta("┴─┘└─┘└┴┘"))
    #
    @pcs "A Modern Toolkit for Analytical Continuation Problems\n" black
    @pcs "Package: " black "$__LIBNAME__\n" magenta
    @pcs "Version: " black "$__VERSION__\n" magenta
    @pcs "Release: " black "$__RELEASE__\n" magenta
    #
    println()
    #
    flush(stdout)
end

"""
    overview()

Print out the overview of ACFlow to the screen.
"""
function overview()
    # Build strings
    str1 = nprocs() == 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # Write the information
    println("[ Overview ]")
    println("Time : ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Para : Using ", nprocs(), str1, str2)
    println("Dirs : ", pwd())
    println("Task : ", query_args())
    println()
    #
    flush(stdout)
end

"""
    goodbye()

Print the goodbye messages to the screen.
"""
function goodbye()
    println("The analytical continuation problem is solved successfully.")
    #
    flush(stdout)
end

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    prompt(msg::String)

Print a stylized ACFlow message to the screen.
"""
function prompt(msg::String)
    print(green("ACFlow > "))
    print(magenta(msg))
    println()
    #
    flush(stdout)
end

#=
### *I/O Operations*
=#

"""
    line_to_array(io::IOStream)

Convert a line (reading from an IOStream) to a string array.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end

#=
### *Math* : *Root Finding*
=#

"""
    secant(func, x0, args...)

It implements the well-known secant algorithm to locate root of a given
polynomial function. Here, `func` means the function, `x0` is the initial
guess, and `args...` denotes the optional parameters for `func`. Please
be careful, `func` must be a single variable function.

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
    newton(fun::Function, guess, kwargs...)

It implements the well-known newton algorithm to locate root of a given
polynomial function. Here, `fun` means the function, `guess` is the initial
solution, and `kwargs...` denotes the optional parameters for `fun`. Please
be careful, `func` is a multiple variable function. It not only returns
the value, but also the jacobian matrix of the function.

See also: [`secant`](@ref).
"""
function newton(fun::Function, guess, kwargs...; maxiter::I64 = 20000, mixing::F64 = 0.5)
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
    trapz(x::AbstractMesh, y::Vector{T})

Perform numerical integration by using the composite trapezoidal rule.

See also: [`simpson`](@ref).
"""
function trapz(x::AbstractMesh, y::Vector{T}) where {T}
    value = dot(x.weight, y)
    return value
end

"""
    trapz(x::Vector{F64}, y::Vector{T}, linear::Bool = false)

Perform numerical integration by using the composite trapezoidal rule.

See also: [`simpson`](@ref).
"""
function trapz(x::Vector{F64}, y::Vector{T}, linear::Bool = false) where {T}
    if linear
        h = x[2] - x[1]
        value = y[1] + y[end] + 2.0 * sum(y[2:end-1])
        value = h * value / 2.0
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
    simpson(x::Vector{F64}, y::Vector{T})

Perform numerical integration by using the simpson rule. Note that the
length of `x` and `y` must be odd numbers. And `x` must be a linear and
uniform mesh.

See also: [`simpson`](@ref).
"""
function simpson(x::Vector{F64}, y::Vector{T}) where {T}
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

abstract type AbstractInterpolation{FT,T} <: AbstractVector{T} end

(interp::AbstractInterpolation)(t::Number) = _interpolate(interp, t)

function munge_data(u::AbstractVector, t::AbstractVector)
    Tu = Base.nonmissingtype(eltype(u))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == length(u)
    non_missing_indices = collect(i for i in 1:length(t) if !ismissing(u[i]) && !ismissing(t[i]))
    newu = Tu.([u[i] for i in non_missing_indices])
    newt = Tt.([t[i] for i in non_missing_indices])
    return newu, newt
end

# Quadratic Interpolation
struct QuadraticInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    t::tType
    function QuadraticInterpolation{FT}(u,t) where {FT}
        new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
    end
end
  
function QuadraticInterpolation(u,t)
    u, t = munge_data(u, t)
    QuadraticInterpolation{true}(u,t)
end

function _interpolate(A::QuadraticInterpolation{<:AbstractVector}, t::Number)
    idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 2))
    i₀, i₁, i₂ = idx, idx + 1, idx + 2
    l₀ = ((t - A.t[i₁])*(t - A.t[i₂]))/((A.t[i₀] - A.t[i₁])*(A.t[i₀] - A.t[i₂]))
    l₁ = ((t - A.t[i₀])*(t - A.t[i₂]))/((A.t[i₁] - A.t[i₀])*(A.t[i₁] - A.t[i₂]))
    l₂ = ((t - A.t[i₀])*(t - A.t[i₁]))/((A.t[i₂] - A.t[i₀])*(A.t[i₂] - A.t[i₁]))
    return A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

# Cubic Spline Interpolation
struct CubicSpline{uType,tType,hType,zType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    t::tType
    h::hType
    z::zType
    function CubicSpline{FT}(u,t,h,z) where {FT}
        new{typeof(u),typeof(t),typeof(h),typeof(z),FT,eltype(u)}(u,t,h,z)
    end
end

function CubicSpline(u,t)
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k+1] - t[k], 1:length(t)-1), 0)
    dl = h[2:n+1]
    d_tmp = 2 .* (h[1:n+1] .+ h[2:n+2])
    du = h[2:n+1]
    tA = Tridiagonal(dl,d_tmp,du)
    d = map(i -> i == 1 || i == n + 1 ? 0 : 6(u[i+1] - u[i]) / h[i+1] - 6(u[i] - u[i-1]) / h[i], 1:n+1)
    z = tA\d
    CubicSpline{true}(u,t,h[1:n+1],z)
end

function _interpolate(A::CubicSpline{<:AbstractVector{<:Number}}, t::Number)
    i = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
    I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
    C = (A.u[i+1]/A.h[i+1] - A.z[i+1]*A.h[i+1]/6)*(t - A.t[i])
    D = (A.u[i]/A.h[i+1] - A.z[i]*A.h[i+1]/6)*(A.t[i+1] - t)
    I + C + D
end

#=
### *Math* : *Einstein Summation Notation*
=#

macro einsum(ex)
    _einsum(ex) # true, false, false
end

function _einsum(expr::Expr, inbounds = true, simd = false, threads = false)
    # Get left hand side (lhs) and right hand side (rhs) of equation
    lhs = expr.args[1]
    rhs = expr.args[2]

    # Get info on the left-hand side
    lhs_arrays, lhs_indices, lhs_axis_exprs = extractindices(lhs)
    length(lhs_arrays) != 1 && throw(ArgumentError(string(
        "Left-hand side of equation contains multiple arguments. Only a single ",
        "referencing expression (e.g. @einsum A[i] = ...) should be used.")))

    # Get info on the right-hand side
    rhs_arrays, rhs_indices, rhs_axis_exprs = extractindices(rhs)
    
    check_index_occurrence(lhs_indices, rhs_indices)

    # remove duplicate indices on left-hand and right-hand side
    # and ensure that the array sizes match along these dimensions
    ###########################################################
    dimension_checks = Expr[]

    # remove duplicate indices on the right hand side
    for i in reverse(eachindex(rhs_indices))
        duplicated = false
        di = rhs_axis_exprs[i]
        
        for j = 1:(i - 1)
            if rhs_indices[j] == rhs_indices[i]
                # found a duplicate
                duplicated = true
                dj = rhs_axis_exprs[j]

                # add dimension check ensuring consistency
                push!(dimension_checks, :(@assert $dj == $di))
            end
        end
        
        for j = eachindex(lhs_indices)
            if lhs_indices[j] == rhs_indices[i]
                dj = lhs_axis_exprs[j]
                if Meta.isexpr(expr, :(:=))
                    # expr.head is :=
                    # infer the size of the lhs array
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

    # remove duplicate indices on the left hand side
    for i in reverse(eachindex(lhs_indices))
        duplicated = false
        di = lhs_axis_exprs[i]

        # don't need to check rhs, already done above

        for j = 1:(i - 1)
            if lhs_indices[j] == lhs_indices[i]
                # found a duplicate
                duplicated = true
                dj = lhs_axis_exprs[j]

                # add dimension check
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
        # infer type of allocated array
        #    e.g. rhs_arrays = [:A, :B]
        #    then the following line produces :(promote_type(eltype(A), eltype(B)))
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

    if threads && !Meta.isexpr(expr, :(=)) && !Meta.isexpr(expr, :(:=))
        throw(ArgumentError(string(
            "Threaded @vielsum can only assign with = or := right now. ",
            "To use ", expr.head, " try @einsum instead.")))
        # could allow :(+=) by simply then removing $lhs = zero($T) line
    end

    if threads && length(lhs_indices) == 0
        # this won't actually cause problems, but won't use threads
        throw(ArgumentError(string(
            "Threaded @vielsum can't assign to a scalar LHS. ",
            "Try @einsum instead.")))
    end

    
    # copy the index expression to modify it; loop_expr is the Expr we'll build loops around
    loop_expr = unquote_offsets!(copy(expr))

    # Nest loops to iterate over the destination variables
    if length(rhs_indices) > 0
        # There are indices on rhs that do not appear in lhs.
        # We sum over these variables.

        if !threads # then use temporaries to write into, as before

            # Innermost expression has form s += rhs
            @gensym s
            loop_expr.args[1] = s
            loop_expr.head = :(+=)

            # Nest loops to iterate over the summed out variables
            loop_expr = nest_loops(loop_expr, rhs_indices, rhs_axis_exprs, simd, false)

            # Prepend with s = 0, and append with assignment
            # to the left hand side of the equation.
            lhs_assignment = Expr(assignment_op, lhs, s)

            loop_expr = quote
                local $s = zero($T)
                $loop_expr
                $lhs_assignment
            end

        else # we are threading, and thus should write directly to lhs array

            loop_expr.args[1] = lhs
            loop_expr.head = :(+=)

            loop_expr = nest_loops(loop_expr, rhs_indices, rhs_axis_exprs, simd, false)

            loop_expr = quote
                $lhs = zero($T)
                $loop_expr
            end
        end

        # Now loop over indices appearing on lhs, if any
        loop_expr = nest_loops(loop_expr, lhs_indices, lhs_axis_exprs, false, threads)
    else
        # We do not sum over any indices, only loop over lhs
        loop_expr.head = assignment_op
        loop_expr = nest_loops(loop_expr, lhs_indices, lhs_axis_exprs, simd, threads)
    end

    if inbounds
        loop_expr = :(@inbounds $loop_expr)
    end

    full_expression = quote
        $type_definition
        $output_definition
        $(dimension_checks...)
        
        # remove let when we drop 0.6 support -- see #31
        let $([lhs_indices; rhs_indices]...)
            $loop_expr
        end

        $(lhs_arrays[1])
    end

    return esc(full_expression)
end

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

"""
    nest_loops(expr, indices, axis_exprs, simd, threads) -> Expr

Construct a nested loop around `expr`, using `indices` in ranges `axis_exprs`.

### Example
```julia-repl
julia> nest_loops(:(A[i] = B[i]), [:i], [:(size(A, 1))], true, false)
quote
    local i
    @simd for i = 1:size(A, 1)
        A[i] = B[i]
    end
end
```
"""
function nest_loops(expr::Expr,
                    index_names::Vector{Symbol}, axis_expressions::Vector{Expr},
                    simd::Bool, threads::Bool)
    isempty(index_names) && return expr

    # Add @simd to the innermost loop, if required
    # and @threads to the outermost loop
    expr = nest_loop(expr, index_names[1], axis_expressions[1],
                     simd, threads && length(index_names) == 1)

    # Add remaining for loops
    for j = 2:length(index_names)
        expr = nest_loop(expr, index_names[j], axis_expressions[j],
                         false, threads && length(index_names) == j)
    end

    return expr
end


"""
    extractindices(expr) -> (array_names, index_names, axis_expressions)

Compute all `index_names` and respective axis calculations of an expression 
involving the arrays with `array_names`. Everything is ordered by first 
occurence in `expr`.

### Examples
```julia-repl
julia> extractindices(:(f(A[i,j,i]) + C[j]))
(Symbol[:A, :C], Symbol[:i, :j, :i, :j], Expr[:(size(A, 1)), :(size(A, 2)), :(size(A, 3)), :(size(C, 1))])
```
"""
extractindices(expr) = extractindices!(expr, Symbol[], Symbol[], Expr[])

function extractindices!(expr::Symbol,
                         array_names::Vector{Symbol},
                         index_names::Vector{Symbol},
                         axis_expressions::Vector{Expr})
    push!(array_names, expr)
    return array_names, index_names, axis_expressions
end

#=
### *Color Tools*
=#

#=
*Remarks* :

The purpose of the following codes is to provide some convenient tools
to output colorful and stylized texts in the terminal. Actually, these
codes are inspried by this repository:

* https://github.com/Aerlinger/AnsiColor.jl

For more information about the ANSI color escape sequences, please check
the following websites further:

* https://stackoverflow.com/questions/4842424/
* https://en.wikipedia.org/wiki/ANSI_escape_code

Note that the macro `@pcs` and functions `prompt()` rely on these codes.
=#

"""
    COLORS

A global dict, which is used to specify the system colors.
"""
const COLORS = Dict{String,I64}(
    "black"          => 0,
    "red"            => 1,
    "green"          => 2,
    "yellow"         => 3,
    "blue"           => 4,
    "magenta"        => 5,
    "cyan"           => 6,
    "white"          => 7,
    "default"        => 9,
    "light_black"    => 60,
    "light_red"      => 61,
    "light_green"    => 62,
    "light_yellow"   => 63,
    "light_blue"     => 64,
    "light_magenta"  => 65,
    "light_cyan"     => 66,
    "light_white"    => 67
)

"""
    MODES

A global dict, which is used to specify the mode for output characters.
"""
const MODES = Dict{String,I64}(
    "default"        => 0,
    "bold"           => 1,
    "underline"      => 4,
    "blink"          => 5,
    "swap"           => 7,
    "hide"           => 8
)

"""
    colorize(c::String, s::String; bg::String = "default", m::String="default")

Return some escape sequences, which will be displayed as colorized texts
in the terminal.
"""
function colorize(c::String, s::String; bg::String = "default", m::String="default")
    C_OFFSET = 30
    B_OFFSET = 40
    "\033[$(MODES[m]);$(C_OFFSET + COLORS[c]);$(B_OFFSET + COLORS[bg])m$(s)\033[0m"
end

"""
    colorize(c::String, s::String; bg::String = "default", m::String="default")

Return some escape sequences, which will be displayed as colorized texts
in the terminal.
"""
function colorize(c::Symbol, s::String; bg::String = "default", m::String="default")
    colorize(string(c), s; bg=bg, m=m)
end

#=
*Remarks* :

The following codes will generate and export dynamically some color
functions, including:

```julia
# For standard colors
black(str::String)
red(str::String)
green(str::String)
yellow(str::String)
blue(str::String)
magenta(str::String)
cyan(str::String)
white(str::String)
```

and their light color versions

```julia
# For light colors
light_black(str::String)
light_red(str::String)
light_green(str::String)
light_yellow(str::String)
light_blue(str::String)
light_magenta(str::String)
light_cyan(str::String)
light_white(str::String)
```

These functions provide some shortcuts to create texts decorated by
special escape sequences. These texts will be show as colorized texts
in the terminal.

### Examples
```julia-repl
julia> println(red("hello world!"))
```
=#

export COLORS
export MODES
export colorize

for k in keys(COLORS)
    f = Symbol(k)
    k == "default" && continue
    @eval ($f)(str::String) = colorize(Symbol($f), str)
    @eval export ($f)
end
