#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

function read_real_data(finput::String, ngrid::I64)
    mesh = zeros(F64, ngrid)
    value = zeros(F64, ngrid)
    error = zeros(F64, ngrid)

    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:3])
            mesh[i] = arr[1]
            value[i] = arr[2]
            error[i] = arr[3]
        end
    end

    return RawData(mesh, value, error)
end

function read_complex_data(finput::String, ngrid::I64)
    mesh = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    open(finput, "r") do fin
        for i = 1:ngrid
            # For test1 and test2
            #arr = parse.(F64, line_to_array(fin)[1:3])
            #mesh[i] = arr[1]
            #value[i] = arr[2] + im * arr[3]
            #error[i] = 0.0001 + 0.0001im

            # For test3 and test4
            #arr = parse.(F64, line_to_array(fin)[1:4])
            #mesh[i] = arr[1]
            #value[i] = 0.0 + im * arr[3]
            #error[i] = arr[4] + im * arr[4]

            # For test5 and test6 and test8
            arr = parse.(F64, line_to_array(fin)[1:4])
            mesh[i] = arr[1]
            value[i] = arr[2] + im * arr[3]
            error[i] = arr[4] + im * arr[4]

            # For test7
            #arr = parse.(F64, line_to_array(fin)[1:3])
            #mesh[i] = arr[1]
            #value[i] = arr[2]
            #error[i] = arr[3]
        end
    end

    return RawData(mesh, value, error)
end

function write_spectrum(am::AbstractMesh, spectra::Vector{F64})
    @assert am.nmesh == length(spectra)
    open("spectra.data", "w") do fout
        for i = 1:am.nmesh
            @printf(fout, "%16.12f %16.12f\n", am[i], spectra[i])
        end
    end
end