#
# Project : Gardenia
# Source  : base.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
#

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})

Solve the analytic continuation problem. The arguments `grid`, `Gval`,
and `Gerr` are the grid, value, and error bar, respectively.

### Arguments
* grid -> Imaginary axis grid for correlators, τ or ωₙ.
* Gval -> Function values for correlators, G(τ) or G(iωₙ).
* Gerr -> Standard deviations for correlators.

### Returns
* mesh -> Real frequency mesh, ω, Vector{F64}.
* Aout -> Spectral function, A(ω), Vector{F64}.
* Gout -> Retarded Green's function, G(ω), Vector{C64}.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T}) where {T}
    return solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, err::T)

Solve the analytic continuation problem. The arguments `grid`, `Gval`,
and `err` are the grid, value, and error bar, respectively.

Here, we just assume that the standard deviations for correlators are
fixed to a constant value, `err`.

### Arguments
* grid -> Imaginary axis grid for correlators, τ or ωₙ.
* Gval -> Function values for correlators, G(τ) or G(iωₙ).
* err  -> Standard deviations for correlators.

### Returns
* mesh -> Real frequency mesh, ω, Vector{F64}.
* Aout -> Spectral function, A(ω), Vector{F64}.
* Gout -> Retarded Green's function, G(ω), Vector{C64}.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, err::T) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, err)
    return solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T})

Solve the analytic continuation problem. The arguments `grid` and `Gval`
are the grid and value, respectively. Furthermore, the error bar is set to
a fixed value `1.0e-4`.

### Arguments
* grid -> Imaginary axis grid for correlators, τ or ωₙ.
* Gval -> Function values for correlators, G(τ) or G(iωₙ).

### Returns
* mesh -> Real frequency mesh, ω, Vector{F64}.
* Aout -> Spectral function, A(ω), Vector{F64}.
* Gout -> Retarded Green's function, G(ω), Vector{C64}.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}) where {T}
    Gerr = similar(Gval)
    err = 1.0e-4
    #
    if T == F64
        fill!(Gerr, err)
    else
        fill!(Gerr, err + err * im)
    end
    #
    return solve(RawData(grid, Gval, Gerr))
end

"""
    solve(rd::RawData)

Solve the analytic continuation problem. The input data are encapsulated
in a `RawData` struct. This function call is the actual interface to the
desired analytic continuation solvers.

### Arguments
* rd -> A RawData struct that contains the grid, correator, and error bar.

### Returns
* mesh -> Real frequency mesh, ω, Vector{F64}.
* Aout -> Spectral function, A(ω) or A(ω) / ω, Vector{F64}.
* Gout -> Retarded Green's function, G(ω), Vector{C64}.

### Examples
```julia
# Setup the configuration file
setup_args("ac.toml")

# Read the parameters
read_param()

# Call the solver.
mesh, Aout, Gout = solve(read_data())
```

See also: [`RawData`](@ref).
"""
function solve(rd::RawData)
    # Return a valid solver object
    function make_solver()
        solver = get_b("solver")

        # Choose suitable solver
        @cswitch solver begin
            @case "MaxEnt"
                return MaxEntSolver()
                break

            @case "BarRat"
                return BarRatSolver()
                break

            @case "NevanAC"
                return NevanACSolver()
                break

            @case "StochAC"
                return StochACSolver()
                break

            @case "StochSK"
                return StochSKSolver()
                break

            @case "StochOM"
                return StochOMSolver()
                break

            @case "StochPX"
                return StochPXSolver()
                break

            @default
                sorry()
                break
        end
    end

    # We just use try...catch block to catch possible exceptions
    # or errors during simulations.
    try
        return solve(make_solver(), rd)
    catch ex
        catch_error()
    end
end

"""
    reprod(am::AbstractMesh, kernel::Matrix{F64}, A::Vector{F64})

Try to reproduce the input data, which can be compared with the raw data
to see whether the analytic continuation is reasonable.

### Arguments
* am -> Real frequency mesh.
* kernel -> The kernel function.
* A -> The calculated spectral function, A(ω) or A(ω) / ω.

### Returns
* G -> Reconstructed correlators, G(τ) or G(iωₙ), Vector{F64}.

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
*Remarks* : *Kramers-Kronig Transformation*

The real and imaginary parts of Green's functions obey the following
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

Try to calculate the real part of the Green's function from its imaginary
part via the Kramers-Kronig relations. The objective of this function is
to get the full retarded Green's function.

### Arguments
* am -> Real frequency mesh.
* A  -> Spectral function at real frequency mesh, A(ω).

### Returns
* G  -> Retarded Green's function, G(ω), Vector{C64}.
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
    setup_param(
        C::Dict{String,Any},
        S::Dict{String,Any},
        reset::Bool = true
    )

Setup the configuration dictionaries via function call. Here `C` contains
parameters for general setup, while `S` contains parameters for selected
analytic continuation solver. If `reset` is true, then the configuration
dictionaries will be reset to their default values at first. Later, `C`
and `S` will be used to customized the dictionaries further.

### Arguments
See above explanations.

### Returns
N/A

See also: [`read_param`](@ref).
"""
function setup_param(
    C::Dict{String,Any},
    S::Dict{String,Any},
    reset::Bool = true
    )
    # _PBASE, _PMaxEnt, _PBarRat, _PNevanAC,
    # _PStochAC, _PStochSK, _PStochOM, and _PStochPX
    #
    # contain the default parameters. If reset is true, they will be used
    # to update the
    #
    # PBASE, PMaxEnt, PBarRat, PNevanAC,
    # PStochAC, PStochSK, PStochOM, and PStochPX
    #
    # dictionaries, respectively.
    reset && begin
        rev_dict_b(_PBASE)
        rev_dict_m(MaxEntSolver(),   _PMaxEnt)
        rev_dict_r(BarRatSolver(),   _PBarRat)
        rev_dict_n(NevanACSolver(), _PNevanAC)
        rev_dict_a(StochACSolver(), _PStochAC)
        rev_dict_k(StochSKSolver(), _PStochSK)
        rev_dict_s(StochOMSolver(), _PStochOM)
        rev_dict_x(StochPXSolver(), _PStochPX)
    end

    rev_dict_b(C)
    #
    solver = get_b("solver")
    #
    @cswitch solver begin
        @case "MaxEnt"
            rev_dict_m(MaxEntSolver(),  S)
            break

        @case "BarRat"
            rev_dict_r(BarRatSolver(),  S)
            break

        @case "NevanAC"
            rev_dict_n(NevanACSolver(), S)
            break

        @case "StochAC"
            rev_dict_a(StochACSolver(), S)
            break

        @case "StochSK"
            rev_dict_k(StochSKSolver(), S)
            break

        @case "StochOM"
            rev_dict_s(StochOMSolver(), S)
            break

        @case "StochPX"
            rev_dict_x(StochPXSolver(), S)
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

### Arguments
N/A

### Returns
N/A

See also: [`setup_param`](@ref).
"""
function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
    chk_dict()
    see_dict()
end

"""
    read_data(only_real_part::Bool = true)

Read data in imaginary axis and return a `RawData` struct. The argument
`only_real_part` is only useful for bosonic cases. If the kernel is
`bsymm` (it means a symmetric bosonic kernel) and the grid is `bfreq` or
`bfrag` (it means Matsubara frequency grid), the function values for
input correators should be real in principle. In these cases, we should
set `only_real_part = true`.

### Arguments
See above explanations.

### Returns
* rd -> Raw input data that is encapsulated in a `RawData` struct.

See also: [`RawData`](@ref).
"""
function read_data(only_real_part::Bool = true)
    finput = get_b("finput")
    ktype = get_b("ktype")
    ngrid = get_b("ngrid")

    function read_dispatcher()
        @cswitch get_b("grid") begin
            @case "ftime"
                return read_real_data(finput, ngrid)
                break

            @case "fpart"
                return read_real_data(finput, ngrid)
                break

            @case "btime"
                return read_real_data(finput, ngrid)
                break

            @case "bpart"
                return read_real_data(finput, ngrid)
                break

            @case "ffreq"
                return read_cmplx_data(finput, ngrid)
                break

            @case "ffrag"
                return read_cmplx_data(finput, ngrid)
                break

            @case "bfreq"
                if ktype == "boson"
                    return read_cmplx_data(finput, ngrid)
                else # ktype == "bsymm"
                    return read_cmplx_data(finput, ngrid, only_real_part)
                end
                break

            @case "bfrag"
                if ktype == "boson"
                    return read_cmplx_data(finput, ngrid)
                else # ktype == "bsymm"
                    return read_cmplx_data(finput, ngrid, only_real_part)
                end
                break

            @default
                sorry()
                break
        end
    end

    # We just use try...catch block to catch possible exceptions
    # or errors during simulations.
    try
        return read_dispatcher()
    catch ex
        catch_error()
    end
end

"""
    make_data(rd::RawData; T::DataType = F64)

Convert `RawData` struct to `GreenData` struct. Note that `RawData` is
provided by the users directly, while `GreenData` is more suitable for
various analytic continuation solvers and algorithms. Note that the
`GreenData` struct is accessed and manipulated by this code internally,
while the `RawData` struct is exposed to the users.

### Arguments
See above explanations.

### Returns
* gd -> A GreenData struct.

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function make_data(rd::RawData; T::DataType = F64)
    ktype = get_b("ktype")
    grid = get_b("grid")
    val = rd.value
    err = rd.error

    if grid == "ffreq" || ( grid == "bfreq" && ktype == "boson" ) ||
       grid == "ffrag" || ( grid == "bfrag" && ktype == "boson" )
        value = T.( vcat(real(val), imag(val)) )
        error = T.( vcat(real(err), imag(err)) )
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    else # grid == "bfreq" && ktype == "bsymm"
         # grid == "bfrag" && ktype == "bsymm"
         # grid == "ftime" || grid == "fpart"
         # grid == "btime" || grid == "bpart"
        value = T.( real(val) )
        error = T.( real(err) )
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end
end

"""
    make_grid(rd::RawData; T::DataType = F64)

Extract grid for input data from a `RawData` struct (`RD`). It will return
a sub-type of the AbstractGrid struct.

### Arguments
See above explanations.

### Returns
* grid -> Imaginary time or imaginary frequency grid.

See also: [`RawData`](@ref), [`AbstractGrid`](@ref).
"""
function make_grid(rd::RawData; T::DataType = F64)
    grid = get_b("grid")
    ngrid = get_b("ngrid")
    β::T = get_b("beta")

    v = T.(rd._grid)
    @assert ngrid == length(v)

    _grid = nothing
    #
    @cswitch grid begin
        @case "ftime"
            _β = v[end]
            @assert abs(_β - β) ≤ 1e-6
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "fpart"
            _grid = FermionicFragmentTimeGrid(β, v)
            break

        @case "btime"
            _β = v[end]
            @assert abs(_β - β) ≤ 1e-6
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "bpart"
            _grid = BosonicFragmentTimeGrid(β, v)
            break

        @case "ffreq"
            _β = 2.0 * π / (v[2] - v[1])
            @assert abs(_β - β) ≤ 1e-6
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
            break

        @case "ffrag"
            _grid = FermionicFragmentMatsubaraGrid(β, v)
            break

        @case "bfreq"
            _β = 2.0 * π / (v[2] - v[1])
            @assert abs(_β - β) ≤ 1e-6
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
            break

        @case "bfrag"
            _grid = BosonicFragmentMatsubaraGrid(β, v)
            break

        @default
            sorry()
            break
    end

    return _grid
end

"""
    make_mesh(; T::DataType = F64)

Try to generate an uniform (linear) or non-uniform (non-linear) mesh for
the spectral function in real axis. Notice that it supports arbitrary
precision mesh. By default, the precision is F64. One can specify the
precision by the argument `T`.

### Arguments
See above explanations.

### Returns
* mesh -> Real frequency mesh. It should be a subtype of AbstractMesh.

See also: [`LinearMesh`](@ref), [`TangentMesh`](@ref), [`LorentzMesh`](@ref).
"""
function make_mesh(; T::DataType = F64)
    # Predefined parameters for mesh generation
    #
    # Note that the parameters `f1` and `cut` are only for the generation
    # of the non-uniform mesh.
    #
    f1::T = 2.1
    cut::T = 0.01

    # Setup parameters according to case.toml
    pmesh = get_b("pmesh")
    #
    if !isa(pmesh, Missing)
        (length(pmesh) == 1) && begin
            Γ, = pmesh
            f1 = Γ
            cut = Γ
        end
    end

    # Get essential parameters
    ktype = get_b("ktype")
    mesh  = get_b("mesh")
    nmesh = get_b("nmesh")
    wmax::T = get_b("wmax")
    wmin::T = get_b("wmin")
    #
    # For bosonic correlators of Hermitian operators, the spectral
    # function is defined in (0, ∞) only.
    if ktype == "bsymm"
        @assert wmin ≥ 0.0
        @assert wmax ≥ 0.0
        @assert wmax > wmin
    end

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

### Arguments
* am -> Real frequency mesh.

### Returns
* model -> A normalized model function.

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
    #
    if !isa(pmodel, Missing)
        (length(pmodel) == 1) && begin Γ, = pmodel end
        (length(pmodel) == 2) && begin Γ, s = pmodel end
        (length(pmodel) == 3) && begin Γ, s₁, s₂ = pmodel end
    end

    # Try to generate the required model
    mtype = get_b("mtype")
    #
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

### Arguments
* am -> Real frequency mesh.
* ag -> Imaginary axis grid.

### Returns
* kernel -> Kernel function, a 2D array, (ntime,nmesh) or (nfreq,nmesh).

See also: [`AbstractMesh`](@ref), [`AbstractGrid`](@ref).
"""
function make_kernel(am::AbstractMesh, ag::AbstractGrid)
    ktype = get_b("ktype")
    grid = get_b("grid")

    @cswitch ktype begin
        @case "fermi"
            @assert grid in ("ftime", "fpart", "ffreq", "ffrag")
            return build_kernel(am, ag)
            break

        @case "boson"
            @assert grid in ("btime", "bpart", "bfreq", "bfrag")
            return build_kernel(am, ag)
            break

        @case "bsymm"
            @assert grid in ("btime", "bpart", "bfreq", "bfrag")
            return build_kernel_symm(am, ag)
            break

        @default
            sorry()
            break
    end
end
