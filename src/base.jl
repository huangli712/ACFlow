#
# Project : Gardenia
# Source  : base.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/01/17
#

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `Gerr` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T}) where {T}
    return solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, err::T)

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `err` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, err::T) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, err)
    return solve(RawData(grid, Gval, Gerr))
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
    return solve(RawData(grid, Gval, Gerr))
end

"""
    solve(rd::RawData)

Solve the analytical continuation problem. The input data are encapsulated
in a `Rawdata` struct.

See also: [`RawData`](@ref).
"""
function solve(rd::RawData)
    solver = get_b("solver")

    # Choose suitable solver
    @cswitch solver begin
        @case "MaxEnt"
            return solve(MaxEntSolver(),  rd)
            break

        @case "StochAC"
            return solve(StochACSolver(), rd)
            break

        @case "StochSK"
            return solve(StochSKSolver(), rd)
            break

        @case "StochOM"
            return solve(StochOMSolver(), rd)
            break

        @case "StochPX"
            return solve(StochPXSolver(), rd)
            break

        @default
            sorry()
            break
    end
end

"""
    reprod(am::AbstractMesh, kernel::Matrix{F64}, A::Vector{F64})

Try to reproduce the input data using the calculated spectrum function
`A`. `kernel` is the kernel function, and `am` is the mesh in which the
spectrum is defined.

See also: [`AbstractMesh`](@ref).
"""
function reprod(am::AbstractMesh, kernel::Matrix{F64}, A::Vector{F64})
    ndim, nmesh = size(kernel)
    @assert nmesh == length(am) == length(A)

    @einsum KA[i,j] := kernel[i,j] * A[j]

    G = zeros(F64, ndim)
    for i = 1:ndim
        G[i] = trapz(am, view(KA, i, :))
    end

    return G
end

#=
*Remarks* : Kramers-Kronig Transformation

The real and imaginary parts of green's functions obey the following
Kramers-Kronig relation.

```math
\begin{equation}
A(\omega) = -\frac{1}{\pi} \mathrm{Im} G(\omega)
\end{equation}
```

```math
\begin{equation}
\mathrm{Re} G(\omega) = \frac{1}{\pi} \mathcal{P}
  \int_{-\infty}^{\infty} d\omega'~
  \frac{\mathrm{Im} G(\omega')}{\omega'-\omega}
\end{equation}
```

```math
\begin{equation}
\mathrm{Re} G(\omega) = \frac{2}{\pi} \mathcal{P}
  \int_{0}^{\infty} d\omega'~
  \frac{\omega'\mathrm{Im} G(\omega')}{\omega'^2 - \omega^2}
\end{equation}
```

So that we can calculate the real parts from the imaginary parts of the
response function, and vice versa.
=#

"""
    kramers(am::AbstractMesh, A::Vector{F64})

Try to calculate the real part of the green's function from its imaginary
part via the Kramers-Kronig relations.
"""
function kramers(am::AbstractMesh, A::Vector{F64})
    nmesh = length(am)
    spectrum = reshape(A, (nmesh,1))
    weight = reshape(am.weight, (nmesh,1))
    w₁ = reshape(am.mesh, (1,nmesh))
    w₂ = reshape(am.mesh, (nmesh,1))

    # For fermionic system
    if am[1] < 0.0
        m = weight .* spectrum ./ (w₁ .- w₂)
    # For bosonic system, the mesh is defined in positive half-axis only.
    else
        m = 2.0 * weight .* w₁ .* spectrum ./ (w₁ .^ 2.0 .- w₂ .^ 2.0)
    end

    # Setup the diagonal elements
    for i = 1:nmesh
        m[i,i] = 0.0
    end

    gre = reshape(sum(m, dims = 1), (nmesh))
    gim = -A * π

    return map((x,y) -> x + im * y, gre, gim)
end

"""
    setup_param(C::Dict{String,Any}, S::Dict{String,Any}, reset::Bool = true)

Setup the configuration dictionaries via function call. Here `C` contains
parameters for general setup, while `S` contains parameters for selected
analytical continuation solver. If `reset` is true, then the configuration
dictionaries will be reset to their default values at first. Later, `C`
`S` will be used to customized the dictionaries further.

See also: [`read_param`](@ref), [`rev_dict`](@ref).
"""
function setup_param(C::Dict{String,Any}, S::Dict{String,Any}, reset::Bool = true)
    # _PBASE, _PMaxEnt, _PStochAC, _PStochSK, _PStochOM, and _PStochPX
    # contain the default parameters. If reset is true, they will be used
    # to update the PBASE, PMaxEnt, PStochAC, PStochSK, PStochOM, and
    # PStochPX dictionaries, respectively.
    reset && begin
        rev_dict(_PBASE)
        rev_dict(MaxEntSolver(),   _PMaxEnt)
        rev_dict(StochACSolver(), _PStochAC)
        rev_dict(StochSKSolver(), _PStochSK)
        rev_dict(StochOMSolver(), _PStochOM)
        rev_dict(StochPXSolver(), _PStochPX)
    end

    rev_dict(C)
    #
    solver = get_b("solver")
    @cswitch solver begin
        @case "MaxEnt"
            rev_dict(MaxEntSolver(),  S)
            break

        @case "StochAC"
            rev_dict(StochACSolver(), S)
            break

        @case "StochSK"
            rev_dict(StochSKSolver(), S)
            break

        @case "StochOM"
            rev_dict(StochOMSolver(), S)
            break

        @case "StochPX"
            rev_dict(StochPXSolver(), S)
            break

        @default
            sorry()
            break
    end
end

"""
    read_param()

Setup the configuration dictionaries via an external file. The valid
format of a configuration file is `toml`.

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
    finput = get_b("finput")
    ktype = get_b("ktype")
    ngrid = get_b("ngrid")

    @cswitch get_b("grid") begin
        @case "ftime"
            return read_real_data(finput, ngrid)
            break

        @case "btime"
            return read_real_data(finput, ngrid)
            break

        @case "ffreq"
            return read_cmplx_data(finput, ngrid)
            break

        @case "bfreq"
            if ktype == "boson"
                return read_cmplx_data(finput, ngrid)
            else
                return read_cmplx_data(finput, ngrid, only_real_part)
            end
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
    ktype = get_b("ktype")
    grid = get_b("grid")
    val = rd.value
    err = rd.error

    if grid == "ffreq" || ( grid == "bfreq" && ktype == "boson" )
        value = vcat(real(val), imag(val))
        error = vcat(real(err), imag(err))
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    else # grid == "bfreq" && ktype == "bsymm"
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
    grid = get_b("grid")
    ngrid = get_b("ngrid")
    v = rd._grid
    @assert ngrid == length(v)

    _grid = nothing
    @cswitch grid begin
        @case "ftime"
            β = v[end]
            @assert abs(β - get_b("beta")) ≤ 1e-6
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "btime"
            β = v[end]
            @assert abs(β - get_b("beta")) ≤ 1e-6
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "ffreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_b("beta")) ≤ 1e-6
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
            break

        @case "bfreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_b("beta")) ≤ 1e-6
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
            break

        @default
            sorry()
            break
    end

    return _grid
end

"""
    make_mesh()

Try to generate an uniform (linear) or non-uniform (non-linear) mesh for
the spectral function in real axis.

See also: [`LinearMesh`](@ref), [`TangentMesh`](@ref), [`LorentzMesh`](@ref).
"""
function make_mesh()
    # Predefined parameters for mesh generation
    #
    # Note that the parameters `f1` and `cut` are only for the generation
    # of the non-uniform mesh.
    #
    f1 = 2.1
    cut = 0.01

    # Setup parameters according to case.toml
    pmesh = get_b("pmesh")
    if !isa(pmesh, Missing)
        (length(pmesh) == 1) && begin
            Γ, = pmesh
            f1 = Γ
            cut = Γ
        end
    end

    # Get essential parameters
    nmesh = get_b("nmesh")
    mesh = get_b("mesh")
    wmax = get_b("wmax")
    wmin = get_b("wmin")

    # Try to generate the required mesh
    @cswitch mesh begin
        @case "linear"
            return LinearMesh(nmesh, wmin, wmax)
            break

        @case "tangent"
            return TangentMesh(nmesh, wmin, wmax, f1)
            break

        @case "lorentz"
            return LorentzMesh(nmesh, wmin, wmax, cut)
            break

        @case "halflorentz"
            return HalfLorentzMesh(nmesh, wmax, cut)
            break

        @default
            sorry()
            break
    end
end

"""
    make_model(am::AbstractMesh)

Try to generate a default model function at given mesh `am` through
various schemes.

See also: [`AbstractMesh`](@ref).
"""
function make_model(am::AbstractMesh)
    # Predefined parameters for model generation
    Γ = 2.0
    s = 2.0
    s₁ = -2.0
    s₂ = +2.0
    fn = "model.inp"

    # Setup parameters according to case.toml
    pmodel = get_b("pmodel")
    if !isa(pmodel, Missing)
        (length(pmodel) == 1) && begin Γ, = pmodel end
        (length(pmodel) == 2) && begin Γ, s = pmodel end
        (length(pmodel) == 3) && begin Γ, s₁, s₂ = pmodel end
    end

    # Try to generate the required model
    mtype = get_b("mtype")
    @cswitch mtype begin
        @case "flat"
            return build_flat_model(am)
            break

        @case "gauss"
            return build_gaussian_model(am, Γ)
            break

        @case "1gauss"
            return build_1gaussian_model(am, Γ, s)
            break

        @case "2gauss"
            return build_2gaussians_model(am, Γ, s₁, s₂)
            break

        @case "lorentz"
            return build_lorentzian_model(am, Γ)
            break

        @case "1lorentz"
            return build_1lorentzian_model(am, Γ, s)
            break

        @case "2lorentz"
            return build_2lorentzians_model(am, Γ, s₁, s₂)
            break

        @case "risedecay"
            return build_risedecay_model(am, Γ)
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
    ktype = get_b("ktype")
    grid = get_b("grid")

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
