#
# Project : Gardenia
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/17
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
 "ac.toml"
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
    newton(fun::Function, guess, kwargs...)

It implements the well-known newton algorithm to locate root of a given
polynomial function. Here, `fun` means the function, `guess` is the initial
solution, and `kwargs...` denotes the required arguments for `fun`. Please
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
    simpson(x::AbstractVector{F64}, y::AbstractVector{T})

Perform numerical integration by using the simpson rule. Note that the
length of `x` and `y` must be odd numbers. And `x` must be a linear and
uniform mesh.

See also: [`simpson`](@ref).
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
directly from https://github.com/PumasAI/DataInterpolations.jl.
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

Preprocess the input data.
"""
function munge_data(u::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    return u, t
end

"""
    munge_data(u::AbstractVector, t::AbstractVector)

Preprocess the input data.
"""
function munge_data(u::AbstractVector, t::AbstractVector)
    Tu = Base.nonmissingtype(eltype(u))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == length(u)
    non_missing_indices = collect(i for i in 1:length(t) if !ismissing(u[i]) && !ismissing(t[i]))
    newu = Tu.([u[i] for i in non_missing_indices])
    newt = Tt.([t[i] for i in non_missing_indices])
    return newu, newt
end

"""
    (interp::AbstractInterpolation)(t::Number)

Support `interp(t)`-like function call, where `interp` is the instance of
any AbstractInterpolation struct.
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
