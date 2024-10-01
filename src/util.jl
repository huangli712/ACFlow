#
# Project : Gardenia
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
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
toolkit is minimizing the dependence on the third-party libraries as
far as possible. Note that the `ACFlow` toolkit relys on the `TOML`
package to parse the *.toml file. Only in v1.6.0 and higher versions,
julia includes the `TOML` package in its standard library.

### Arguments
N/A

### Returns
N/A
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

### Arguments
* x -> Filename of configuration file.

### Returns
* ARGS -> Global variable.

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

### Arguments
N/A

### Returns
* x -> ARGS[1], where ARGS is a global variable.

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
### *Error Handler*
=#

"""
    trace_error(io)

Write exceptions or errors to terminal or external file.

### Arguments
* io -> Output stream.

### Returns
N/A

See also: [`catch_error`](@ref).
"""
function trace_error(io)
    # current_exceptions() will return the stack of exceptions
    # currently being handled.
    for (exc, btrace) in current_exceptions()
        Base.showerror(io, exc, btrace)
        println(io)
    end
end

"""
    catch_error()

Catch the thrown exceptions or errors, print them to the terminal or
external file (`err.out`).

### Arguments
N/A

### Returns
N/A

### Examples
```julia
try
    return solve(read_data())
catch ex
    catch_error()
end
```

See also: [`trace_error`](@ref).
"""
function catch_error()
    # For REPL case, error messages are written to terminal
    if isinteractive()
        println(red("ERROR: "), magenta("The stacktrace is shown below"))
        trace_error(stdout)
    # For standard case, error messages will be written into err.out.
    else
        println("ERROR: The stacktrace is saved in err.out")
        open("err.out", "a") do fio
            trace_error(fio)
        end
    end
end

#=
### *Colorful Outputs*
=#

"""
    welcome()

Print out the welcome messages to the screen.

### Arguments
N/A

### Returns
N/A
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

### Arguments
N/A

### Returns
N/A
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

### Arguments
N/A

### Returns
N/A
"""
function goodbye()
    println("The analytic continuation problem is solved successfully.")
    println("Current Time : ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    #
    flush(stdout)
end

"""
    sorry()

Print an error message to the screen.

### Arguments
N/A

### Returns
N/A
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    prompt(msg::String)

Print a stylized ACFlow message to the screen.

### Arguments
* msg -> Message that need to be printed.

### Returns
N/A
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

### Arguments
* io -> An IOStream struct.

### Returns
* arr -> An array  of String.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.

### Arguments
* str -> A String.

### Returns
* ass -> An array of String.

### Examples
```julia-repl
julia> str = "Hello World!"
"Hello World!"

julia> line_to_array(str)
2-element Vector{SubString{String}}:
 "Hello"
 "World!"
```
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
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
    colorize(
        c::String,
        s::String;
        bg::String = "default",
        m::String = "default"
    )

Return some escape sequences, which will be displayed as colorized texts
in the terminal.

### Arguments
* c  -> Color names.
* s  -> The string that want to be printed.
* bg -> Background color.
* m  -> Output mode.

### Returns
* See above explanations.
"""
function colorize(
    c::String,
    s::String;
    bg::String = "default",
    m::String = "default"
    )
    C_OFFSET = 30
    B_OFFSET = 40
    "\033[$(MODES[m]);$(C_OFFSET + COLORS[c]);$(B_OFFSET + COLORS[bg])m$(s)\033[0m"
end

"""
    colorize(
        c::Symbol,
        s::String;
        bg::String = "default",
        m::String = "default"
    )

Return some escape sequences, which will be displayed as colorized texts
in the terminal.

### Arguments
* c  -> Color names.
* s  -> The string that want to be printed.
* bg -> Background color.
* m  -> Output mode.

### Returns
* See above explanations.
"""
function colorize(
    c::Symbol,
    s::String;
    bg::String = "default",
    m::String = "default"
    )
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
