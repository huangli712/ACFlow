#
# Project : Gardenia
# Source  : rfa.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/07/26
#

#=
### *Customized Structs* : *BarRat Solver*
=#

mutable struct  BarRatContext
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

function init(S::BarRatSolver, rd::RawData)
end

function run(brc::BarRatContext)
end

function last(brc::BarRatContext)
end
