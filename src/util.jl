#
# Project : Gardenia
# Source  : util.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/20
#

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
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

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

function gamma(x::F64)
    ccall(("tgamma", libm), F64, (F64,), x)
end
