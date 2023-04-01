#!/usr/bin/env julia

#
# This script is used to generate non-uniform mesh according to the
# given spectral function. It will launch only 1 process.
#
# Usage:
#
#     $ gmesh.jl ac.toml
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Printf
using ACFlow

function read_spectrum(fn::String = "Aout.data")
    nmesh = get_b("nmesh")

    mesh = zeros(F64, nmesh)
    image = zeros(F64, nmesh)

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

    return mesh, image
end

function write_fmesh(fn::String = "fmesh.inp", fmesh::Vector{F64})
    open(fn, "w") do fout
        for i in eachindex(fmesh)
            @printf(fout, "%8i %16.12f\n", i, fmesh[i])
        end
    end
end

function generate_mesh(mesh::Vector{F64}, image::Vector{F64})
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

generate_mesh(read_spectrum()...)

goodbye()
