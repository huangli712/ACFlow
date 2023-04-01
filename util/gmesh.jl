#!/usr/bin/env julia

#
# This script is used to generate a non-uniform very dense mesh according
# to the given spectral function. It will launch only 1 process.
#
# Usage:
#
#     $ gmesh.jl ac.toml
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Printf
using ACFlow

"""
    read_spectrum(fn::String = "Aout.data")

Try to read the spectral function from the present directory. It will
return the spectrum and the corresponding mesh in real axis.
"""
function read_spectrum(fn::String = "Aout.data")
    # Get essential parameter from the configuration file
    nmesh = get_b("nmesh")

    # Allocate memories for mesh and spectrum
    mesh = zeros(F64, nmesh)
    image = zeros(F64, nmesh)

    # Extract the spectral data
    if isfile(fn)
        open(fn, "r") do fin
            for i = 1:nmesh
                arr = line_to_array(fin)
                mesh[i], image[i] = parse.(F64, arr)
            end
        end
    else
        error("Sorry, $fn does not exist!")
    end

    # Return the spectral data
    return mesh, image
end

"""
    write_fmesh(fmesh::Vector{F64}, fn::String = "fmesh.inp")

Write the generated mesh to the file `fmesh.inp`.
"""
function write_fmesh(fmesh::Vector{F64}, fn::String = "fmesh.inp")
    open(fn, "w") do fout
        for i in eachindex(fmesh)
            @printf(fout, "%8i %16.12f\n", i, fmesh[i])
        end
    end
end

"""
    generate_fmesh(mesh::Vector{F64}, image::Vector{F64})


"""
function generate_fmesh(mesh::Vector{F64}, image::Vector{F64})
    nfine = 100000
    Asum = cumsum(abs.(image))
    f = LinearInterpolation(mesh, Asum)

    vmax = maximum(Asum)
    vmin = minimum(Asum)
    Amesh = collect(LinRange(vmin, vmax, nfine))
    fmesh = f.(Amesh)
    return fmesh
end

welcome()
overview()
read_param()

fmesh = generate_fmesh(read_spectrum()...)
write_fmesh(fmesh, "test.inp")

goodbye()
