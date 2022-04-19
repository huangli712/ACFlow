#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/04/19
#

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `Gerr` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T}) where {T}
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, err::T)

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `err` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector, Gval::Vector{T}, err::T) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, err)
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T})

Solve the analytical continuation problem. The arguments `grid` and `Gval`
are the grid and value, respectively. Furthermore, the error bar is set to
a fixed value `1.0e-4`.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, 1.0e-4)
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(rd::RawData)

Solve the analytical continuation problem. The input data are encapsulated
in a `Rawdata` struct.

See also: [`RawData`](@ref).
"""
function solve(rd::RawData)
    solver = get_c("solver")

    @cswitch solver begin
        @case "MaxEnt"
            solve(MaxEntSolver(), rd)
            break

        @case "StochAC"
            solve(StochACSolver(), rd)
            break

        @case "StochOM"
            solve(StochOMSolver(), rd)
            break

        @default
            sorry()
            break
    end
end

"""
    reprod(kernel::Matrix{F64}, am::AbstractMesh, A::Vector{F64})

Try to reproduce the input data using the calculated spectrum function
`A`. `kernel` is the kernel function, and `am` is the mesh in which the
spectrum is defined.

See also: [`AbstractMesh`](@ref).
"""
function reprod(kernel::Matrix{F64}, am::AbstractMesh, A::Vector{F64})
    ndim, nmesh = size(kernel)
    @assert nmesh == length(am) == length(A)

    @einsum K[i,j] := kernel[i,j] * A[j]

    G = zeros(F64, ndim)
    for i = 1:ndim
        G[i] = trapz(am, view(K, i, :))
    end

    return G
end

#=
*Remarks* : Kramers-Kronig transformation
```math
\begin{equation}
A(\omega) = -\frac{1}{\pi} \mathrm{Im} G(\omega)
\end{equation}
```

```math
\begin{equation}
\mathrm{Re} G(\omega) = \frac{1}{\pi} \mathcal{P} 
\int_{-\infty}^\infty d\omega~\frac{\mathrm{Im} G(\omega')}{\omega'-\omega}
\end{equation}
```
=#

"""
    kramers(am::AbstractMesh, A::Vector{F64})

The Kramers-Kronig relations provide a way to get the real part from the
imaginary part.
"""
function kramers(am::AbstractMesh, A::Vector{F64})
    nmesh = length(am)
    spectrum = reshape(A, (nmesh,1))
    weight = reshape(am.weight, (nmesh,1))
    w₁ = reshape(am.mesh, (1,nmesh))
    w₂ = reshape(am.mesh, (nmesh,1))

    m = weight .* spectrum ./ (w₁ .- w₂)
    for i = 1:nmesh
        m[i,i] = 0.0
    end

    gre = reshape(sum(m, dims = 1), (nmesh))
    gim = -A * π

    return gre, gim
end

"""
    setup_param(C::Dict{String,Any}, S::Dict{String,Any})

Setup the configuration dictionaries via function call. Here `C` contains
parameters for general setup, while `S` contains parameters for selected
analytical continuation solver.

See also: [`read_param`](@ref), [`rev_dict`](@ref).
"""
function setup_param(C::Dict{String,Any}, S::Dict{String,Any})
    rev_dict(C)

    solver = get_c("solver")
    @cswitch solver begin
        @case "MaxEnt"
            rev_dict(MaxEntSolver(), S)
            break

        @case "StochAC"
            rev_dict(StochACSolver(), S)
            break

        @case "StochOM"
            rev_dict(StochOMSolver(), S)
            break

        @default
            sorry()
            break
    end
end

"""
    read_param()

Setup the configuration dictionaries via a external file. The valid format
of a configuration file is `toml`.

See also: [`read_param`](@ref).
"""
function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
    chk_dict()
end

"""
    read_data(only_real_part::Bool = true)

Read data in imaginary axis and return a `RawData` struct.

See also: [`RawData`](@ref).
"""
function read_data(only_real_part::Bool = true)
    finput = get_c("finput")
    ngrid = get_c("ngrid")

    @cswitch get_c("grid") begin
        @case "ftime"
            return read_real_data(finput, ngrid)
            break

        @case "btime"
            return read_real_data(finput, ngrid)
            break

        @case "ffreq"
            return read_complex_data(finput, ngrid)
            break

        @case "bfreq"
            return read_complex_data(finput, ngrid, only_real_part)
            break

        @default
            sorry()
            break
    end
end

"""
    make_data(rd::RawData)

Convert `RawData` struct to `GreenData` struct. Note that `RawData` is
provided by the users directly, while `GreenData` is more suitable for
various analytical continuation solvers and algorithms. Note that the
`GreenData` struct is accessed and manipulated by this code internally,
while the `RawData` struct is exposed to the users.

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function make_data(rd::RawData)
    grid = get_c("grid")
    val = rd.value
    err = rd.error

    if grid == "ffreq"
        value = vcat(real(val), imag(val))
        error = vcat(real(err), imag(err))
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    else
        value = real(val)
        error = real(err)
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end
end

"""
    make_grid(rd::RawData)

Extract grid for input data from a `RawData` struct. It will return a
sub-type of the AbstractGrid struct.

See also: [`RawData`](@ref), [`AbstractGrid`](@ref).
"""
function make_grid(rd::RawData)
    grid = get_c("grid")
    ngrid = get_c("ngrid")
    v = rd._grid
    @assert ngrid == length(v)

    _grid = nothing
    @cswitch grid begin
        @case "ftime"
            β = v[end]
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "btime"
            β = v[end]
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "ffreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
            break

        @case "bfreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
            break

        @default
            sorry()
            break
    end

    return _grid
end

"""
    make_mesh(f1::F64 = 2.1)

Try to generate an uniform (linear) or non-uniform (non-linear) mesh for
the calculated spectrum. Note that the argument `f1` is only for the
generation of the non-uniform mesh.

See also: [`LinearMesh`](@ref), [`TangentMesh`](@ref).
"""
function make_mesh(f1::F64 = 2.1)
    nmesh = get_c("nmesh")
    mesh = get_c("mesh")
    wmax = get_c("wmax")
    wmin = get_c("wmin")

    @cswitch mesh begin
        @case "linear"
            return LinearMesh(nmesh, wmin, wmax)
            break

        @case "tangent"
            return TangentMesh(nmesh, wmin, wmax, f1)
            break

        @default
            sorry()
            break
    end
end

"""
    make_model(am::AbstractMesh, Γ::F64 = 2.0, fn::String = "model.data")

Try to generate default model function through various schemes. Note that
the argument `Γ` is for the gauss-like model function, while the argument
`fn` is for the user-supplied model function.

See also: [`AbstractMesh`]
"""
function make_model(am::AbstractMesh, Γ::F64 = 2.0, fn::String = "model.data")
    mtype = get_c("mtype")

    @cswitch mtype begin
        @case "flat"
            return build_flat_model(am)
            break

        @case "gauss"
            return build_gaussian_model(am, Γ)
            break

        @case "file"
            return build_file_model(am, fn)
            break

        @case "func"
            sorry()
            break

        @default
            sorry()
            break
    end
end

"""
    make_kernel(am::AbstractMesh, ag::AbstractGrid)

Try to generate various kernel functions.

See also: [`AbstractMesh`](@ref), [`AbstractGrid`](@ref).
"""
function make_kernel(am::AbstractMesh, ag::AbstractGrid)
    ktype = get_c("ktype")
    grid = get_c("grid")

    @cswitch ktype begin
        @case "fermi"
            return build_kernel(am, ag)
            break

        @case "boson"
            return build_kernel(am, ag)
            break

        @case "bsymm"
            @assert grid in ("btime", "bfreq")
            return build_kernel_symm(am, ag)
            break

        @default
            sorry()
            break
    end
end
