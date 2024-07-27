#
# Project : Gardenia
# Source  : rfa.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/07/27
#

#=
### *Customized Structs* : *BarRat Solver*
=#

"""
    Barycentric (type)

Barycentric representation of a rational function.

# Fields
- `node`: the nodes of the rational function
- `value`: the values of the rational function
- `weight`: the weights of the rational function
- `wf`: the weighted values of the rational function
"""
struct Barycentric{T} <: Function
    nodes::Vector{C64}
    values::Vector{C64}
    weights::Vector{C64}
    w_times_f::Vector{C64}
    function Barycentric{T}(
        node::AbstractVector{C64},
        value::AbstractVector{C64},
        weight::AbstractVector{C64},
        wf::AbstractVector{C64} = value.*weight
        ) where {T <: AbstractFloat}
        @assert length(node) == length(value) == length(weight) == length(wf)
        new{T}(node, value, weight, wf)
    end
end

"""
    Barycentric(node, value, weight, wf=value.*weight)

Construct a `Barycentric` rational function.

# Arguments
- `node::AbstractVector`: interpolation nodes
- `value::AbstractVector`: values at the interpolation nodes
- `weight::AbstractVector`: barycentric weights
- `wf::AbstractVector`: weights times values (optional)

# Examples
```jldoctest
julia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])
Barycentric function with 3 nodes and values:
    1.0=>1.0,  2.0=>2.0,  3.0=>3.0

julia> r(1.5)
1.5
```
"""
function Barycentric(
    node::Vector{S}, value::Vector{S}, weight::Vector{S}, wf=value.*weight
    ) where {S<:C64}
    return Barycentric{F64}(node, value, weight, wf)
end

"""
    BarRatContext

Mutable struct. It is used within the BarRat solver only.

### Members

* Gᵥ   -> Input data for correlator.
* grid -> Grid for input data.
* mesh -> Mesh for output spectrum.
"""
mutable struct  BarRatContext
    Gᵥ   :: Vector{C64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    ℬ    :: Union{Missing,Barycentric{F64}}
end

#=
### *Global Drivers*
=#

"""
    solve(S::BarRatSolver, rd::RawData)

Solve the analytic continuation problem by the Barycentric rational
function method.
"""
function solve(S::BarRatSolver, rd::RawData)
    println("[ BarRat ]")
    #
    brc = init(S, rd)
    run(brc)
    Aout, Gout = last(brc)
    #
    return brc.mesh.mesh, Aout, Gout
end

"""
    init(S::BarRatSolver, rd::RawData)

Initialize the BarRat solver and return a BarRatContext struct.
"""
function init(S::BarRatSolver, rd::RawData)
    # Prepera input data
    Gᵥ = rd.value
    println("Postprocess input data: ", length(Gᵥ), " points")

    # Prepare grid for input data
    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    return BarRatContext(Gᵥ, grid, mesh, missing)
end

function run(brc::BarRatContext)
    r = aaa(brc.grid.ω * im, brc.Gᵥ)
    @show r
    @show poles(r)
    @show typeof(r)
    brc.ℬ = r
end

function last(brc::BarRatContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Calculate full response function on real axis and write them
    _G = brc.ℬ.(brc.mesh.mesh)
    fwrite && write_complete(brc.mesh, _G)

    # Calculate and write the spectral function
    Aout = -imag.(_G) ./ π
    fwrite && write_spectrum(brc.mesh, Aout)

    # Regenerate the input data and write them
    kernel = make_kernel(brc.mesh, brc.grid)
    G = reprod(brc.mesh, kernel, Aout)
    fwrite && write_backward(brc.grid, G)

    return Aout, _G
end

"""
    poles(r)

Return the poles of the rational function `r`.
"""
function poles(r::Barycentric{T}) where T
    w = weights(r)
    nonzero = @. !iszero(w)
    z, w = nodes(r)[nonzero], w[nonzero]
    m = length(w)
    B = diagm( [zero(T); ones(T, m)] )
    E = [zero(T) transpose(w); ones(T, m) diagm(z) ];
    pol = []  # put it into scope
    try
        pol = filter( isfinite, eigvals(E, B) )
    catch
        # generalized eigen not available in extended precision, so:
        λ = filter( z->abs(z)>1e-13, eigvals(E\B) )
        pol = 1 ./ λ
    end
    return pol
end

"weights(r) returns the weights of the rational interpolant `r` as a vector."
weights(r::Barycentric) = r.weights

"nodes(r) returns the nodes of the rational interpolant `r` as a vector."
nodes(r::Barycentric) = r.nodes

"""
    r(z)
    evaluate(r, z)

Evaluate the rational function at `z`.
"""

(r::Barycentric)(z) = evaluate(r, z)
function evaluate(r::Barycentric, z::Number)
    if isinf(z)
        return sum(r.w_times_f) / sum(r.weights)
    end
    k = findfirst(z .== r.nodes)
    if isnothing(k)         # not at a node
        C = @. 1 / (z - r.nodes)
        return sum(C .* r.w_times_f) / sum(C .* r.weights)
    else                    # interpolation at node
        return r.values[k]
    end
end

"""
    aaa(z, y)

Adaptively compute a rational interpolant.

# Arguments

## discrete mode
- `z::AbstractVector{<:Number}`: interpolation nodes
- `y::AbstractVector{<:Number}`: values at nodes

## continuous mode
- `f::Function`: function to approximate on the interval [-1,1]

# Keyword arguments
- `max_degree::Integer=150`: maximum numerator/denominator degree to use
- `float_type::Type=Float64`: floating point type to use for the computation
- `tol::Real=1000*eps(float_type)`: tolerance for stopping
- `lookahead::Integer=10`: number of iterations to determines stagnation
- `stats::Bool=false`: return convergence statistics

# Returns
- `r::Barycentric`: the rational interpolant
- `stats::NamedTuple`: convergence statistics, if keyword `stats=true`

# Examples
```julia-repl
julia> z = 1im * range(-10, 10, 500);

julia> y = @. exp(z);

julia> r = aaa(z, y);

julia> degree(r)   # both numerator and denominator
12

julia> first(nodes(r), 4)
4-element Vector{ComplexF64}:
 0.0 - 6.272545090180361im
 0.0 + 9.43887775551102im
 0.0 - 1.1022044088176353im
 0.0 + 4.909819639278557im

julia> r(1im * π / 2)
-2.637151617496356e-15 + 1.0000000000000002im
```

See also [`approximate`](@ref) for approximating a function on a curve or region.
"""
function aaa(z::AbstractVector{<:Number}, y::AbstractVector{<:Number};
    max_degree = 150, float_type = Float64, tol = 1000*eps(float_type),
    lookahead = 10, stats = false
    )

    @assert float_type <: AbstractFloat
    T = float_type
    fmax = norm(y, Inf)    # for scaling
    m = length(z)
    iteration = NamedTuple[]
    err = T[]
    besterr, bestidx, best = Inf, NaN, nothing

    # Allocate space for Cauchy matrix, Loewner matrix, and residual
    C = similar(z, (m, m))
    L = similar(z, (m, m))
    R = complex(zeros(size(z)))

    ȳ = sum(y) / m
    s, idx = findmax(abs(y - ȳ) for y in y)
    push!(err, s)

    # The ordering of nodes matters, while the order of test points does not.
    node_index = Int[]
    push!(node_index, idx)
    test_index = Set(1:m)
    delete!(test_index, idx)

    n = 0 # number of poles
    while true
        n += 1
        σ = view(z, node_index)
        fσ = view(y, node_index)
        # Fill in matrices for the latest node
        @inbounds @fastmath for i in test_index
            δ = z[i] - σ[n]
            # δ can be zero if there are repeats in z
            C[i, n] = iszero(δ) ? 1 / eps() : 1 / δ
            L[i, n] = (y[i] - fσ[n]) * C[i, n]
        end
        istest = collect(test_index)
        _, _, V = svd( view(L, istest, 1:n) )
        w = V[:, end]    # barycentric weights

        CC = view(C, istest, 1:n)
        num = CC * (w.*fσ)
        den = CC * w
        @. R[istest] = y[istest] - num / den
        push!(err, norm(R, Inf))
        push!(iteration, (; weights=w, active=copy(node_index)))

        if (Base.last(err) < besterr)
            besterr, bestidx, best = Base.last(err), length(iteration), Base.last(iteration)
        end

        # Are we done?
        if (besterr <= tol*fmax) ||
            (n == max_degree + 1) ||
            ((length(iteration) - bestidx >= lookahead) && (besterr < 1e-2*fmax))
            break
        end

        _, j = findmax(abs, R)
        push!(node_index, j)
        delete!(test_index, j)
        R[j] = 0
    end

    idx, w = best.active, best.weights
    r = Barycentric(z[idx], y[idx], w)
    if stats
        return r, (;err, iteration)
    else
        return r
    end
end
