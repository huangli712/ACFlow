#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/19
#

#=
### *Input Functions*
=#

"""
    read_real_data(finput::AbstractString, ngrid::I64)

Read input data. This function is used for imaginary-time data. The input
file should contain three columns. The first column is the imaginary-time
grid, the second column is the value, the third column is the standard
deviation σ. Here, `ngrid` specifies the number of grid points.

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
    @assert nrows == ngrid
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
    @assert nrows == ngrid
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
    read_cmplx_data(finput::AbstractString, ngrid::I64, only_real_part::Bool)

Read input data. This function is used for Matsubara frequency data. The
input file only contains three columns. The first column is the Matsubara
frequency grid, the second column is the real part or imaginary part of
the data (which is specified by the argument `only_real_part`), and the
third column is the standard deviation σ. This function is for bosonic
correlation function.

See also: [`read_real_data`](@ref).
"""
function read_cmplx_data(finput::AbstractString, ngrid::I64, only_real_part::Bool)
    # Allocate memory
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    # We have to determine the number of columns and rows at first.
    dlm = readdlm(finput)
    nrows, ncols = size(dlm)
    @assert nrows == ngrid
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

Write spectrum A(ω) to `Aout.data`. The grid is defined in `am`, and
the spectrum is contained in `Aout`.
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
axis, `G` is the calculated green's function data. Note that its real
part is obtained via the so-called Kramers-Kronig transformation.

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
judge whether the obtained optimal α parameter is reasonable.

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

Write the default model function to `model.data`.
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
    write_hamiltonian(α_vec::Vector{F64}, Uα::Vector{F64})

Write `α-U(α)` data to `hamil.data`, which could be used to judge whether
the obtained optimal α parameter is reasonable. This function is only
useful for the `StochAC` solver.
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
    write_pole(Pᵥ::Vector{Vector{I64}},
               Aᵥ::Vector{Vector{F64}},
               χ²::Vector{F64},
               fmesh::AbstractMesh)

Write poles' positions and amplitudes to `pole.data`. This function is
only useful for the `StochPX` solver.
"""
function write_pole(Pᵥ::Vector{Vector{I64}},
                    Aᵥ::Vector{Vector{F64}},
                    χ²::Vector{F64},
                    fmesh::AbstractMesh)
    ntry = length(Pᵥ)

    open("pole.data", "w") do fout
        for i = 1:ntry
            println(fout, "# Try: ", i, "  χ²: ", χ²[i])
            for j in eachindex(Pᵥ[i])
                p = Pᵥ[i][j]
                @printf(fout, "%4i %8i %16.12f %16.12f\n", j, p, fmesh[p], Aᵥ[i][j])
            end
            println(fout)
            println(fout)
        end
    end
end

"""
    write_probability(α_vec::Vector{F64}, p_vec::Vector{F64})

Write `p(α)` data to `prob.data`. This function is only useful for the
`MaxEnt` solver (`bryan` algorithm).
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
