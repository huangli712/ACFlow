#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/07
#

"""
    read_real_data
"""
function read_real_data(finput::String, ngrid::I64)
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
    read_complex_data
"""
function read_complex_data(finput::String, ngrid::I64; ncols::I64 = 4)
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
    read_complex_data
"""
function read_complex_data(finput::String, ngrid::I64, only_real_part::Bool)
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
    write_spectrum
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
    write_chi2
"""
function write_chi2(α_vec::Vector{F64}, χ_vec::Vector{F64})
    @assert length(α_vec) == length(χ_vec)
    open("chi2.data", "w") do fout
        α_vec = log10.(α_vec)
        χ_vec = log10.(χ_vec)
        for i = 1:length(α_vec)
            @printf(fout, "%16.12f %16.12f\n", α_vec[i], χ_vec[i])
        end
    end
end

"""
    write_probability
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
    write_reprod
"""
function write_reprod(ag::AbstractGrid, G::Vector{F64})
    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng || ngrid * 2 == ng
    
    if ngrid == ng
        open("Grep.data", "w") do fout
            for i = 1:ngrid
                @printf(fout, "%16.12f %16.12f\n", ag[i], G[i])
            end
        end
    else
        open("Grep.data", "w") do fout
            for i = 1:ngrid
                @printf(fout, "%16.12f %16.12f %16.12f\n", ag[i], G[i], G[i+ngrid])
            end
        end
    end
end
