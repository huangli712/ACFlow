#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/11
#

"""
    read_real_data(finput::AbstractString, ngrid::I64)

Read input data. This function is used for imaginary-time data. The input
file should contain three columns. The first column is the imaginary-time
grid, the second column is the value, the third column is the standard
deviation σ. Here, `ngrid` specifies the number of grid points.

See also: [`read_complex_data`](@ref).
"""
function read_real_data(finput::AbstractString, ngrid::I64)
    _grid = zeros(F64, ngrid)
    value = zeros(F64, ngrid)
    error = zeros(F64, ngrid)

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
    read_complex_data(finput::AbstractString, ngrid::I64; ncols::I64 = 4)

Read input data. This function is used for Matsubara frequency data. The
input should contain four columns or five columns. The first column is
the Matsubara freqency grid, the second and third columns are the values
(real part and imaginary part), the four and fifth columns are the standard
deviations for the real and imaginary parts, respectively. If there are
only four columns, it means that the real and imaginary parts share the
same standard deviations.

See also: [`read_time_data`](@ref).
"""
function read_complex_data(finput::AbstractString, ngrid::I64; ncols::I64 = 4)
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    @assert ncols in (4, 5)

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
    read_complex_data(finput::AbstractString, ngrid::I64, only_real_part::Bool)

Read input data. This function is used for Matsubara frequency data. The
input file only contains three columns. The first column is the Matsubara
frequency grid, the second column is the real part or imaginary part of
the data (which is specified by the argument `only_real_part`), and the
third column is the standard deviation.

See also: [`read_time_data`](@ref).
"""
function read_complex_data(finput::AbstractString, ngrid::I64, only_real_part::Bool)
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

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

"""
    write_spectrum(am::AbstractMesh, Aout::Vector{F64})

Write spectrum A(ω) to `Aout.data`. The grid is contained in `am`, and
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

Write α-resolved spectrum A(ω) to `Aout.data.alpha`. The grid is contained
in `am`, the α-resolved spectrum is contained in `Aout`, `αₗ` is the list
for the α parameters.
"""
function write_spectrum(am::AbstractMesh, αₗ::Vector{F64}, Aout::Array{F64,2})
    nmesh, nalph = size(Aout)
    @assert nmesh == length(am)
    @assert nalph == length(αₗ)

    
end

"""
    write_chi2(α_vec::Vector{F64}, χ²_vec::Vector{F64})

Write `log10(α)-log10(χ²)` data to `chi2.data`.
"""
function write_chi2(α_vec::Vector{F64}, χ²_vec::Vector{F64})
    @assert length(α_vec) == length(χ²_vec)
    open("chi2.data", "w") do fout
        α_vec = log10.(α_vec)
        χ²_vec = log10.(χ²_vec)
        for i = 1:length(α_vec)
            @printf(fout, "%16.12f %16.12f\n", α_vec[i], χ²_vec[i])
        end
    end
end

"""
    write_probability(α_vec::Vector{F64}, p_vec::Vector{F64})

Write `p(α)` data to `prob.data`. This function is only useful for MaxEnt
bryan algorithm.
"""
function write_probability(α_vec::Vector{F64}, p_vec::Vector{F64})
    @assert length(α_vec) == length(p_vec)
    open("prob.data", "w") do fout
        p_vec = -p_vec ./ trapz(α_vec, p_vec)
        α_vec = log10.(α_vec)
        for i = 1:length(α_vec)
            @printf(fout, "%16.12f %16.12f\n", α_vec[i], p_vec[i])
        end
    end
end

"""
    write_reprod(ag::AbstractGrid, G::Vector{F64})

We can use the calculated spectrum in real axis to reproduce the input
data in imaginary axis. This function will write the reproduced data to
`repr.data`, which can be compared with the original data.
"""
function write_reprod(ag::AbstractGrid, G::Vector{F64})
    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng || ngrid * 2 == ng

    if ngrid == ng
        open("repr.data", "w") do fout
            for i = 1:ngrid
                @printf(fout, "%16.12f %16.12f\n", ag[i], G[i])
            end
        end
    else
        open("repr.data", "w") do fout
            for i = 1:ngrid
                @printf(fout, "%16.12f %16.12f %16.12f\n", ag[i], G[i], G[i+ngrid])
            end
        end
    end
end
