#
# Project : Gardenia
# Source  : spx.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/03
#

#=
### *Customized Structs* : *StochPX Solver*
=#

mutable struct StochPXElement
    P :: Vector{I64}
    A :: Vector{F64}
end

mutable struct StochPXContext
    Gᵥ    :: Vector{F64}
    Gᵧ    :: Vector{F64}
    σ¹    :: Vector{F64}
    allow :: Vector{I64}
    grid  :: AbstractGrid
    mesh  :: AbstractMesh
    Aout  :: Vector{F64}
end

"""
    solve(S::StochPXSolver, rd::RawData)

Solve the analytical continuation problem by the stochastic
pole expansion.
"""
function solve(S::StochPXSolver, rd::RawData)
    nmesh = get_b("nmesh")

    println("[ StochPX ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochPX)
        #
        # Launch the task
        𝐹 = Future[]
        for i = 1:nworkers()
            𝑓 = @spawnat i + 1 prun(S, p1, p2, MC, SE, SC)
            push!(𝐹, 𝑓)
        end
        #
        # Wait and collect the solutions
        sol = []
        for i = 1:nworkers()
            wait(𝐹[i])
            push!(sol, fetch(𝐹[i]))
        end
        #
        # Average the solutions
        Aout = zeros(F64, nmesh)
        for i in eachindex(sol)
            a = sol[i]
            @. Aout = Aout + a / nworkers()
        end
        #
        # Postprocess the solutions
        Gout = last(SC, Aout)

    # Sequential version
    else
        Aout = run(MC, SE, SC)
        Gout = last(SC, Aout)

    end

    return SC.mesh.mesh, Aout, Gout
end

function init(S::StochPXSolver, rd::RawData)
    # Initialize possible constraints. The allow array contains all the
    # possible indices for poles.
    allow = constraints(S)

    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(S, MC.rng, allow)
    println("Randomize Monte Carlo configurations")

end

function run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
end

function prun(S::StochPXSolver,
            p1::Dict{String,Vector{Any}},
            p2::Dict{String,Vector{Any}},
            MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
end

function average(step::F64, SC::StochPXContext)
end

function last(SC::StochPXContext, Aout::Array{F64,2})
end

function warmup()
end

function sample()
end

function measure()
end

#=
### *Service Functions*
=#

"""
    init_mc(S::StochPXSolver)

Try to create a StochPXMC struct.

See also: [`StochPXMC`](@ref).
"""
function init_mc(S::StochPXSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Pacc = 0
    Ptry = 0
    Aacc = 0
    Atry = 0

    MC = StochPXMC(rng, Pacc, Ptry, Aacc, Atry)

    return MC
end

function init_element(S::StochPXSolver)
end

function init_iodata(S::StochPXSolver, rd::RawData)
end

"""
    constraints(S::StochPXSolver)

Try to implement the constrained stochastic pole expansion. This
function will return a collection. It contains all the allowable indices.

See also: [`StochPXSolver`](@ref).
"""
function constraints(S::StochPXSolver)
    exclude = get_b("exclude")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_x("nfine")

    allow = I64[]
    mesh = collect(LinRange(wmin, wmax, nfine))

    # Go through the fine linear mesh and check each mesh point.
    # Is is excluded ?
    for i in eachindex(mesh)
        is_excluded = false
        #
        if !isa(exclude, Missing)
            for j in eachindex(exclude)
                if exclude[j][1] ≤ mesh[i] ≤ exclude[j][2]
                    is_excluded = true
                    break
                end
            end
        end
        #
        if !is_excluded
            push!(allow, i)
        end
    end

    return allow
end

function try_move_p()
end

function try_move_a()
end
