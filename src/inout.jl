#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/30
#

#=
### *Input Functions*
=#

"""
    read_real_data(finput::AbstractString, ngrid::I64)

Read input data. This function is used for imaginary time data. The input
file should contain three columns. The first column is the imaginary time
grid, the second column is the value, the third column is the standard
deviation σ. Here, `ngrid` specifies the number of grid points.

### Arguments
* finput -> Filename for the input data.
* ngrid  -> Number of grid points.

### Returns
* rd -> A RawData struct.

See also: [`read_cmplx_data`](@ref).
"""
function read_real_data(finput::AbstractString, ngrid::I64)
    # Allocate memory
    _grid = zeros(F64, ngrid)
    value = zeros(F64, ngrid)
    error = zeros(F64, ngrid)

    # We have to determine the number of columns and rows at first.
    dlm = readdlm(finput)
    nrows, ncols = size(dlm)
    @assert nrows ≥ ngrid
    @assert ncols == 3

    # Read and parse the data
    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:3])
            _grid[i] = arr[1]
            value[i] = arr[2]
            error[i] = arr[3]
        end
    end

    return RawData(_grid, value, error)
end

"""
    read_cmplx_data(finput::AbstractString, ngrid::I64)

Read input data. This function is used for Matsubara frequency data. The
input should contain four columns or five columns. The first column is
the Matsubara freqency grid, the second and third columns are the values
(real part and imaginary part), the four and fifth columns are the standard
deviations σ for the real and imaginary parts, respectively. If there are
only four columns, it means that the real and imaginary parts share the
same standard deviations.

### Arguments
* finput -> Filename for the input data.
* ngrid  -> Number of grid points.

### Returns
* rd -> A RawData struct.

See also: [`read_real_data`](@ref).
"""
function read_cmplx_data(finput::AbstractString, ngrid::I64)
    # Allocate memory
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    # We have to determine the number of columns and rows at first.
    dlm = readdlm(finput)
    nrows, ncols = size(dlm)
    @assert nrows ≥ ngrid
    @assert ncols in (4, 5)

    # Read and parse the data
    open(finput, "r") do fin
        for i = 1:ngrid
            if ncols == 4
                arr = parse.(F64, line_to_array(fin)[1:4])
                _grid[i] = arr[1]
                value[i] = arr[2] + im * arr[3]
                error[i] = arr[4] + im * arr[4]
            else
                arr = parse.(F64, line_to_array(fin)[1:5])
                _grid[i] = arr[1]
                value[i] = arr[2] + im * arr[3]
                error[i] = arr[4] + im * arr[5]
            end
        end
    end

    return RawData(_grid, value, error)
end

"""
    read_cmplx_data(
        finput::AbstractString,
        ngrid::I64,
        only_real_part::Bool
    )

Read input data. This function is used for Matsubara frequency data. The
input file only contains three columns. The first column is the Matsubara
frequency grid, the second column is the real part or imaginary part of
the data (which is specified by the argument `only_real_part`), and the
third column is the standard deviation σ. This function is for bosonic
correlation function.

### Arguments
* finput -> Filename for the input data.
* ngrid  -> Number of grid points.
* only_real_part -> See above explanations.

### Returns
* rd -> A RawData struct.

See also: [`read_real_data`](@ref).
"""
function read_cmplx_data(
    finput::AbstractString,
    ngrid::I64,
    only_real_part::Bool
    )
    # Allocate memory
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    # We have to determine the number of columns and rows at first.
    dlm = readdlm(finput)
    nrows, ncols = size(dlm)
    @assert nrows ≥ ngrid
    @assert ncols == 3

    # Read and parse the data
    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:3])
            _grid[i] = arr[1]
            if only_real_part
                value[i] = arr[2]
                error[i] = arr[3]
            else
                value[i] = im * arr[2]
                error[i] = im * arr[3]
            end
        end
    end

    return RawData(_grid, value, error)
end

#=
### *Output Functions*
=#

"""
    write_spectrum(am::AbstractMesh, Aout::Vector{F64})

Write spectrum A(ω) or A(ω) / ω to `Aout.data`. The grid is defined in
`am`, and the spectral data are contained in `Aout`.

Note that for the MaxEnt, StochAC, StochSK, and StochOM solvers, `Aout`
is actually A(ω) / ω, instead of A(ω). However, for the BarRat, NevanAC,
and StochPX solvers, `Aout` is just A(ω).

### Arguments
* am   -> Real frequency mesh.
* Aout -> Spectral function. See above explanations.

### Returns
N/A
"""
function write_spectrum(am::AbstractMesh, Aout::Vector{F64})
    @assert length(am) == length(Aout)

    open("Aout.data", "w") do fout
        for i in eachindex(am)
            @printf(fout, "%16.12f %16.12f\n", am[i], Aout[i])
        end
    end
end

"""
    write_spectrum(am::AbstractMesh, αₗ::Vector{F64}, Aout::Array{F64,2})

Write α-resolved spectrum A(ω) to `Aout.data.alpha`. The grid is defined
in `am`, the α-resolved spectrum is contained in `Aout`, `αₗ` is the list
for the α parameters. This function is called by the `StochAC` solver.

Note that `Aout` is actually A(ω) / ω, instead of A(ω).

### Arguments
* am   -> Real frequency mesh.
* αₗ   -> List for α parameters.
* Aout -> α-dependent spectral function.

### Returns
N/A
"""
function write_spectrum(am::AbstractMesh, αₗ::Vector{F64}, Aout::Array{F64,2})
    nmesh, nalph = size(Aout)
    @assert nmesh == length(am)
    @assert nalph == length(αₗ)

    for i in eachindex(αₗ)
        open("Aout.data.alpha_$i", "w") do fout
            println(fout, "# $i : α = ", αₗ[i])
            for j in eachindex(am)
                @printf(fout, "%16.12f %16.12f\n", am[j], Aout[j,i])
            end
        end
    end
end

"""
    write_backward(ag::AbstractGrid, G::Vector{F64})

We can use the calculated spectrum in real axis to reproduce the input
data in imaginary axis. This function will write the reproduced data to
`repr.data`, which can be compared with the original data. Here, `G` is
the reproduced data.

### Arguments
* ag -> Grid for input data, τ or iωₙ.
* G  -> Reconstructed Green's function, G(τ) or G(iωₙ).

### Returns
N/A

See also: [`reprod`](@ref).
"""
function write_backward(ag::AbstractGrid, G::Vector{F64})
    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng || ngrid * 2 == ng

    # The reproduced data are defined in imaginary time axis.
    if ngrid == ng
        open("repr.data", "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f %16.12f\n", ag[i], G[i])
            end
        end
    # The reproduced data are defined in Matsubara frequency axis.
    else
        open("repr.data", "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f %16.12f %16.12f\n", ag[i], G[i], G[i+ngrid])
            end
        end
    end
end

"""
    write_complete(am::AbstractMesh, G::Vector{C64})

Write the full data at real axis to `Gout.data`. `am` denotes the real
axis, `G` is the calculated Green's function data. Note that its real
part is obtained via the so-called Kramers-Kronig transformation.

### Arguments
* am -> Real frequency mesh, ω.
* G  -> Retarded Green's function, G(ω).

### Returns
N/A

See also: [`kramers`](@ref).
"""
function write_complete(am::AbstractMesh, G::Vector{C64})
    @assert length(am) == length(G)

    open("Gout.data", "w") do fout
        for i in eachindex(am)
            z = G[i]
            @printf(fout, "%16.12f %16.12f %16.12f\n", am[i], real(z), imag(z))
        end
    end
end

"""
    write_misfit(α_vec::Vector{F64}, χ²_vec::Vector{F64})

Write `log10(α)-log10(χ²)` data to `chi2.data`, which could be used to
judge whether the obtained optimal α parameter is reasonable. It is used
by the `MaxEnt` solver only.

### Arguments
* α_vec  -> List for α parameters.
* χ²_vec -> α-dependent goodness-of-fit functional.

### Returns
N/A

See also: [`write_goodness`](@ref).
"""
function write_misfit(α_vec::Vector{F64}, χ²_vec::Vector{F64})
    @assert length(α_vec) == length(χ²_vec)

    open("chi2.data", "w") do fout
        _α = log10.(α_vec)
        _χ² = log10.(χ²_vec)
        for i in eachindex(α_vec)
            @printf(fout, "%16.12f %16.12f\n", _α[i], _χ²[i])
        end
    end
end

"""
    write_goodness(Θ_vec::Vector{F64}, χ²_vec::Vector{F64})

Write `log10(Θ)-log10(χ²)` data to `goodness.data`, which could be used
to judge whether the obtained optimal Θ parameter is reasonable. This
function is only useful for the `StochSK` solver.

### Arguments
* Θ_vec  -> List for Θ parameters.
* χ²_vec -> Θ-dependent goodness-of-fit functional.

### Returns
N/A

See also: [`write_misfit`](@ref).
"""
function write_goodness(Θ_vec::Vector{F64}, χ²_vec::Vector{F64})
    @assert length(Θ_vec) == length(χ²_vec)

    open("goodness.data", "w") do fout
        _Θ = log10.(Θ_vec)
        _χ² = log10.(χ²_vec)
        for i in eachindex(Θ_vec)
            if !isinf(_Θ[i]) && !isinf(_χ²[i])
                @printf(fout, "%16.12f %16.12f\n", _Θ[i], _χ²[i])
            end
        end
    end
end

"""
    write_model(am::AbstractMesh, D::Vector{F64})

Write the default model function to `model.data`. This function is usually
for the `MaxEnt` solver.

### Arguments
* am -> Real frequency mesh, ω.
* D  -> Default model, m(ω).

### Returns
N/A
"""
function write_model(am::AbstractMesh, D::Vector{F64})
    @assert length(am) == length(D)

    open("model.data", "w") do fout
        for i in eachindex(am)
            @printf(fout, "%16.12f %16.12f\n", am[i], D[i])
        end
    end
end

"""
    write_prony(𝑁ₚ::I64, Γₚ::Vector{C64}, Ωₚ::Vector{C64})

Write Prony approximation to the input correlator. This information can
be used to reconstruct or interpolate the correlator. This function is
only useful for the `BarRat` solver.

### Arguments
* 𝑁ₚ -> Number of nodes for Prony approximation.
* Γₚ -> Nodes for Prony approximation, ``γ_i``.
* Ωₚ -> Weights for Prony approximation, ``w_i``.

### Returns
N/A
"""
function write_prony(𝑁ₚ::I64, Γₚ::Vector{C64}, Ωₚ::Vector{C64})
    open("prony.data", "w") do fout
        println(fout, "# Prony Approximation")
        #
        println(fout, "# 𝑁ₚ :")
        @printf(fout, "%4i\n", 𝑁ₚ)
        #
        println(fout, "# Γₚ :")
        for i in eachindex(Γₚ)
            z = Γₚ[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
        #
        println(fout, "# Ωₚ :")
        for i in eachindex(Ωₚ)
            z = Ωₚ[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
    end
end

"""
    write_prony(ag::AbstractGrid, G::Vector{C64})

Write Matsubara Green's function approximated by the Prony approximation.
This function is only useful for the `BarRat` solver.

### Arguments
* ag -> Grid for input data, iωₙ.
* G  -> Approximated Green's function, G(iωₙ).

### Returns
N/A
"""
function write_prony(ag::AbstractGrid, G::Vector{C64})
    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng

    open("Gprony.data", "w") do fout
        for i in eachindex(ag)
            z = G[i]
            @printf(fout, "%16.12f %16.12f %16.12f\n", ag[i], real(z), imag(z))
        end
    end
end

"""
    write_barycentric(
        nodes::Vector{C64},
        values::Vector{C64},
        weights::Vector{C64}
    )

Write barycentric rational function approximation to the input correlator.
This information can be used to reconstruct or interpolate the correlator.
This function is only useful for the `BarRat` solver.

### Arguments
* nodes   -> Nodes of the rational function, ``z_i``.
* values  -> Values of the rational function, ``r(z_i)``.
* weights -> Weights of the rational function, ``w_i``.

### Returns
N/A
"""
function write_barycentric(
    nodes   :: Vector{C64},
    values  :: Vector{C64},
    weights :: Vector{C64}
    )
    w_times_f = values .* weights

    open("barycentric.data", "w") do fout
        println(fout, "# Barycentric Rational Function Approximation")
        #
        println(fout, "# nodes :")
        for i in eachindex(nodes)
            z = nodes[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
        #
        println(fout, "# values :")
        for i in eachindex(values)
            z = values[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
        #
        println(fout, "# weights :")
        for i in eachindex(weights)
            z = weights[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
        #
        println(fout, "# w_times_f :")
        for i in eachindex(w_times_f)
            z = w_times_f[i]
            @printf(fout, "%4i %16.12f %16.12f\n", i, real(z), imag(z))
        end
    end
end

"""
    write_hamiltonian(α_vec::Vector{F64}, Uα::Vector{F64})

Write `α-U(α)` data to `hamil.data`, which could be used to judge whether
the obtained optimal α parameter is reasonable. This function is only
useful for the `StochAC` solver.

### Arguments
* α_vec -> List for α parameters.
* Uα    -> α-dependent Hamiltonian.

### Returns
N/A
"""
function write_hamiltonian(α_vec::Vector{F64}, Uα::Vector{F64})
    @assert length(α_vec) == length(Uα)

    open("hamil.data", "w") do fout
        for i in eachindex(α_vec)
            @printf(fout, "%16.8f %16.12f\n", α_vec[i], Uα[i])
        end
    end
end

"""
    write_passed(passed::Vector{I64}, med::F64, αgood::F64)

Write indices of selected solutions which should be used to calculate
the averaged spectrum. Here, `passed` means the indices, `med` is the
median value of χ², and `αgood` is the factor that is used to filter
the solutions. This function is only useful for the `StochOM` and the
`StochPX` solvers.

### Arguments
* passed -> Indices for selected solutions.
* med    -> Median value of χ².
* αgood  -> Predefined parameter used to filter the solutions (spectra).

### Returns
N/A
"""
function write_passed(passed::Vector{I64}, med::F64, αgood::F64)
    open("passed.data", "w") do fout
        npass = length(passed)
        println(fout, "# Count: ", npass, "  Median: ", med, "  αgood: ", αgood)
        for i in eachindex(passed)
            @printf(fout, "%4i %8i\n", i, passed[i])
        end
    end
end

"""
    write_pole(
        Pᵥ::Vector{Vector{I64}},
        Aᵥ::Vector{Vector{F64}},
        𝕊ᵥ::Vector{Vector{F64}},
        χ²ᵥ::Vector{F64},
        fmesh::AbstractMesh
    )

Write positions, amplitudes, and signs of poles to `pole.data`. This
function is only useful for the `StochPX` solver.

Note that the `util/ppole.jl` script can read poles' information from
the `pole.data` file, and construct new spectral functions with different
`eta` parameters.

### Arguments
* Pᵥ    -> Positions of the poles.
* Aᵥ    -> Amplitudes of the poles.
* 𝕊ᵥ    -> Signs of the poles.
* χ²ᵥ   -> Goodness-of-fit functionals for all the solutions.
* fmesh -> A dense mesh for the poles.

### Returns
N/A
"""
function write_pole(
    Pᵥ::Vector{Vector{I64}},
    Aᵥ::Vector{Vector{F64}},
    𝕊ᵥ::Vector{Vector{F64}},
    χ²ᵥ::Vector{F64},
    fmesh::AbstractMesh
    )
    ntry = length(Pᵥ)

    open("pole.data", "w") do fout
        for i = 1:ntry
            println(fout, "# Try: ", i, "  χ²: ", χ²ᵥ[i])
            #
            for j in eachindex(Pᵥ[i])
                ℙ = Pᵥ[i][j]
                𝔸 = Aᵥ[i][j]
                𝕊 = 𝕊ᵥ[i][j]
                ω = fmesh[ℙ]
                @printf(fout, "%4i %8i %16.12f %16.12f %6.2f\n", j, ℙ, ω, 𝔸, 𝕊)
            end
            #
            println(fout)
            println(fout)
        end
    end
end

"""
    write_probability(α_vec::Vector{F64}, p_vec::Vector{F64})

Write `p(α)` data to `prob.data`. This function is only useful for the
`MaxEnt` solver (`bryan` algorithm).

### Arguments
* α_vec -> List for α parameters.
* p_vec -> α-dependent probabilities.

### Returns
N/A
"""
function write_probability(α_vec::Vector{F64}, p_vec::Vector{F64})
    @assert length(α_vec) == length(p_vec)

    open("prob.data", "w") do fout
        _p = -p_vec ./ trapz(α_vec, p_vec)
        _α = log10.(α_vec)
        for i in eachindex(α_vec)
            @printf(fout, "%16.12f %16.12f\n", _α[i], _p[i])
        end
    end
end

"""
    write_statistics(MC::StochACMC)

Write Monte Carlo statistical information for the `StochAC` solver. Note
that the `StochAC` solver is based on a stochastic approach.

### Arguments
* MC -> A StochACMC struct.

### Returns
N/A

See also: [`PStochAC`](@ref), [`StochACMC`](@ref).
"""
function write_statistics(MC::StochACMC)
    nalph = get_a("nalph")

    open("stat.data", "w") do fout
        println(fout, "# Move statistics:")
        for i = 1:nalph
            @printf(fout, "α %3i -> %16.12f\n", i, MC.Macc[i] / MC.Mtry[i])
        end
        println(fout)
        println(fout, "# Swap statistics:")
        for i = 1:nalph
            @printf(fout, "α %3i -> %16.12f\n", i, MC.Sacc[i] / MC.Stry[i])
        end
    end
end

"""
    write_statistics(MC::StochSKMC)

Write Monte Carlo statistical information for the `StochSK` solver. Note
that the `StochSK` solver is based on a stochastic approach.

### Arguments
* MC -> A StochSKMC struct.

### Returns
N/A

See also: [`PStochSK`](@ref), [`StochSKMC`](@ref).
"""
function write_statistics(MC::StochSKMC)
    open("stat.data", "w") do fout
        println(fout, "# Move S statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Sacc)
        @printf(fout, "try    -> %16i    \n", MC.Stry)
        @printf(fout, "prob   -> %16.12f \n", MC.Sacc / MC.Stry)
        println(fout)
        println(fout, "# Move P statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Pacc)
        @printf(fout, "try    -> %16i    \n", MC.Ptry)
        @printf(fout, "prob   -> %16.12f \n", MC.Pacc / MC.Ptry)
        println(fout)
        println(fout, "# Move Q statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Qacc)
        @printf(fout, "try    -> %16i    \n", MC.Qtry)
        @printf(fout, "prob   -> %16.12f \n", MC.Qacc / MC.Qtry)
    end
end

"""
    write_statistics(MC::StochOMMC)

Write Monte Carlo statistical information for the `StochOM` solver. Note
that the `StochOM` solver is based on a stochastic approach.

### Arguments
* MC -> A StochOMMC struct.

### Returns
N/A

See also: [`PStochOM`](@ref), [`StochOMMC`](@ref).
"""
function write_statistics(MC::StochOMMC)
    open("stat.data", "w") do fout
        println(fout, "# Update statistics:")
        @printf(fout, "insert box -> %16.12f\n", MC.Macc[1] / MC.Mtry[1])
        @printf(fout, "remove box -> %16.12f\n", MC.Macc[2] / MC.Mtry[2])
        @printf(fout, "shift  box -> %16.12f\n", MC.Macc[3] / MC.Mtry[3])
        @printf(fout, "width  box -> %16.12f\n", MC.Macc[4] / MC.Mtry[4])
        @printf(fout, "height box -> %16.12f\n", MC.Macc[5] / MC.Mtry[5])
        @printf(fout, "split  box -> %16.12f\n", MC.Macc[6] / MC.Mtry[6])
        @printf(fout, "merge  box -> %16.12f\n", MC.Macc[7] / MC.Mtry[7])
    end
end

"""
    write_statistics(MC::StochPXMC)

Write Monte Carlo statistical information for the `StochPX` solver. Note
that the `StochPX` solver is based on a stochastic approach.

### Arguments
* MC -> A StochPXMC struct.

### Returns
N/A

See also: [`PStochPX`](@ref), [`StochPXMC`](@ref).
"""
function write_statistics(MC::StochPXMC)
    open("stat.data", "w") do fout
        println(fout, "# Move S statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Sacc)
        @printf(fout, "try    -> %16i    \n", MC.Stry)
        @printf(fout, "prob   -> %16.12f \n", MC.Sacc / MC.Stry)
        println(fout)
        println(fout, "# Move P statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Pacc)
        @printf(fout, "try    -> %16i    \n", MC.Ptry)
        @printf(fout, "prob   -> %16.12f \n", MC.Pacc / MC.Ptry)
        println(fout)
        println(fout, "# Move A statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Aacc)
        @printf(fout, "try    -> %16i    \n", MC.Atry)
        @printf(fout, "prob   -> %16.12f \n", MC.Aacc / MC.Atry)
        println(fout)
        println(fout, "# Move X statistics:")
        @printf(fout, "accept -> %16i    \n", MC.Xacc)
        @printf(fout, "try    -> %16i    \n", MC.Xtry)
        @printf(fout, "prob   -> %16.12f \n", MC.Xacc / MC.Xtry)
    end
end
