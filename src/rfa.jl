#
# Project : Gardenia
# Source  : rfa.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/07/27
#

#=
### *Customized Structs* : *BarRat Solver*
=#

"""
    BarRatContext

Mutable struct. It is used within the BarRat solver only.

### Members

* Gᵥ   -> Input data for correlator.
* grid -> Grid for input data.
* mesh -> Mesh for output spectrum.
"""
mutable struct  BarRatContext
    Gᵥ   :: Vector{C64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
end

#=
### *Global Drivers*
=#

"""
    solve(S::BarRatSolver, rd::RawData)

Solve the analytic continuation problem by the Barycentric rational
function method.
"""
function solve(S::BarRatSolver, rd::RawData)
    println("[ BarRat ]")
    #
    brc = init(S, rd)
    run(brc)
    Aout, Gout = last(brc)
    #
    return brc.mesh.mesh, Aout, Gout
end

"""
    init(S::BarRatSolver, rd::RawData)

Initialize the BarRat solver and return a BarRatContext struct.
"""
function init(S::BarRatSolver, rd::RawData)
    # Prepera input data
    Gᵥ = rd.value
    println("Postprocess input data: ", length(Gᵥ), " points")

    # Prepare grid for input data
    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    return BarRatContext(Gᵥ, grid, mesh)
end

function run(brc::BarRatContext)
end

function last(brc::BarRatContext)
end
