#!/usr/bin/env julia

#
# This script is used to generate a non-uniform and very dense mesh
# according to the given spectral function `Aout.data`. The mesh will
# be stored in `fmesh.inp`. The stochastic analytic continuation and
# the stochastic pole expansion can utilize this mesh to resolve the
# fine structure of the spectrum automatically. This script will launch
# only 1 process.
#
# Usage:
#
#     $ gmesh.jl ac.toml
#

using Printf
using ACFlow

"""
    read_spectrum(fn::String = "Aout.data")

Try to read the spectral function from the present directory. It will
return the spectrum and the corresponding mesh in real axis.

### Arguments
* fn -> Filename for the spectral function.

### Returns
* mesh -> Real freqency mesh, ω.
* Aout -> Spectral function, abs(A(ω)).
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
    #
    # For bosonic systems and matrix-valued Green's functions, the spectra
    # could exhibit negative weights. We need the absolute values.
    return mesh, abs.(image)
end

"""
    write_fmesh(fmesh::Vector{F64}, fn::String = "fmesh.inp")

Write the generated mesh to the file `fmesh.inp`. This file can be used by
some stochastic analytic continuation methods, such as the StochPX solver,
to perform constrained sampling.

### Arguments
* fmesh -> Dynamical mesh.
* fn -> Filename for the output mesh.

### Returns
N/A
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

### Arguments
* mesh -> Real freqency mesh, ω.
* image -> Spectral function, abs(A(ω)).

### Returns
* fmesh -> Generated mesh, which should be written into `fmesh.inp`.
"""
function generate_fmesh(mesh::Vector{F64}, image::Vector{F64})
    # Get the analytic continuation solver
    solver = get_b("solver")

    # Get number of points from the configuration file
    nfine = 10000
    #
    @cswitch solver begin
        @case "MaxEnt"
            sorry()
            break

        @case "BarRat"
            sorry()
            break

        @case "NevanAC"
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

    # Create the linear interpolator, 𝑦 = mesh and 𝑥 = Asum.
    𝑓 = LinearInterpolation(mesh, Asum)

    # Interpolate it. Now Amesh is the new 𝑥, and fmesh is the new 𝑦.
    fmesh = 𝑓.(Amesh)

    # Return the generated mesh
    return fmesh
end

welcome()
overview()
read_param()
write_fmesh( generate_fmesh( read_spectrum()... ) )
goodbye()
