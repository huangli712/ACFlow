#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/06
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
    write_chidata
"""
function write_chidata(am::AbstractMesh)
end
