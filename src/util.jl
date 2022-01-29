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

@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
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

#function myfun(a, b)
#    #return a ^ 2.0 + b * a + 1.0
#    return a ^ 3.0 - b
#end

#s = secant(myfun, 1.0, 8.0)
#println(s)

function secant(func, x0, args)
    eps = 1.0e-4
    maxiter = 50
    tol = 1.48e-8
    funcalls = 0
    p0 = 1.0 * x0
    p1 = x0 * (1.0 + eps)
    if p1 ≥ 0.0
        p1 = p1 + eps
    else
        p1 = p1 - eps
    end

    q0 = func(p0, args)
    funcalls = funcalls + 1
    q1 = func(p1, args)
    funcalls = funcalls + 1

    if abs(q1) < abs(q0)
        p0, p1 = p1, p0
        q0, q1 = q1, q0
    end

    for itr = 1:maxiter
        if q1 == q0
            if p1 != p0
                error("tolerance is reached!")
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
        q1 = func(p1, args)
        funcalls = funcalls + 1
    end
end

function newton(fun::Function, guess, kwargs...)
    max_iter = 20000
    mixing = 0.5
    counter = 0
    result = nothing

    function apply(prop::Vector{T}, f::Vector{T}, J::Matrix{T}) where {T}
        resid = nothing
        step = 1.0
        limit = 1e-4
    
        try
            resid = - pinv(J) * f
        catch
            resid = zeros(F64, length(prop))
        end
    
        if any(x -> x > limit, abs.(prop))
            ratio = abs.(resid ./ prop)
            max_ratio = maximum( ratio[ abs.(prop) .> limit ] )
            if max_ratio > 1.0
                step = 1.0 / max_ratio
            end
        end
    
        result = prop + step .* resid
    
        return result
    end

    props = []
    reals = []

    f, J = fun(guess, kwargs...)
    init = apply(guess, f, J)
    push!(props, guess)
    push!(reals, init)

    while true
        counter = counter + 1

        prop = props[end] + mixing * (reals[end] - props[end])
        f, J = fun(prop, kwargs...)
        result = apply(prop, f, J)

        push!(props, prop)
        push!(reals, result)
        
        any(isnan.(result)) && error("Got NaN!")

        if counter > max_iter || maximum( abs.(result - prop) ) < 1.e-4
            break
        end
    end

    counter > max_iter && error("The max_iter is too small!")

    return result, counter
end