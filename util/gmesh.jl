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
    return mesh, abs.(image)
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

Try to generate a very dense and non-uniform mesh according to the given
mesh and spectral function.
"""
function generate_fmesh(mesh::Vector{F64}, image::Vector{F64})
    # Get the analytical continuation solver
    solver = get_b("solver")

    # Get number of points from the configuration file
    nfine = 10000
    @cswitch solver begin
        @case "MaxEnt"
            sorry()
            break

        @case "StochAC"
            nfine = get_a("nfine")
            break

        @case "StochSK"
            nfine = get_k("nfine")
            break

        @case "StochOM"
            sorry()
            break

        @case "StochPX"
            nfine = get_x("nfine")
            break

        @default
            sorry()
            break
    end

    # Accumulate the spectrum, evaluate its maximum and minimum, and then
    # use them to create a very dense linear grid.
    Asum = cumsum(abs.(image))
    vmax = maximum(Asum)
    vmin = minimum(Asum)
    Amesh = collect(LinRange(vmin, vmax, nfine))

    # Create the linear interpolator, y = mesh and x = Asum.
    f = LinearInterpolation(mesh, Asum)

    # Interpolate it. Now Amesh is the new x, and fmesh is the new y.
    fmesh = f.(Amesh)

    # Return the generated mesh
    return fmesh
end

welcome()
overview()
read_param()

write_fmesh( generate_fmesh( read_spectrum()... ) )

goodbye()
