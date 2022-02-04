#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

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

function read_complex_data(finput::String, ngrid::I64)
    _grid = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:4])
            _grid[i] = arr[1]
            value[i] = arr[2] + im * arr[3]
            error[i] = arr[4] + im * arr[4]
        end
    end

    return RawData(_grid, value, error)
end

function write_spectrum(am::AbstractMesh, spectra::Vector{F64})
    @assert am.nmesh == length(spectra)
    open("spectra.data", "w") do fout
        for i = 1:am.nmesh
            @printf(fout, "%16.12f %16.12f\n", am[i], spectra[i])
        end
    end
end