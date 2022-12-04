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
    fmesh :: AbstractMesh
    χ²    :: Vector{F64}
    Pᵥ    :: Vector{Vector{I64}}
    Aᵥ    :: Vector{Vector{F64}}
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
        Aout, Gout, Gr = run(MC, SE, SC)
        last(SC, Aout, Gout, Gr)

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

    Gᵥ, σ¹ = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    fmesh = calc_fmesh(S)
    Gᵧ = deepcopy(Gᵥ)

    χ², Pᵥ, Aᵥ = init_context(S)
    #@show typeof(χ²), typeof(Pᵥ), typeof(Aᵥ)
    #error()
    SC = StochPXContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, fmesh, χ², Pᵥ, Aᵥ)

    return MC, SE, SC
end

function run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")

    for t = 1:ntry
        println("Try: ", t)

        reset_mc(MC)
        reset_element(MC.rng, SC.allow, SE)
        reset_context(t, SE, SC)
    
        for i = 1:nstep
            sample(t, MC, SE, SC)
            if SC.χ²[t] < 1e-6
                println("try :", t, " ", SC.χ²[t])
                break
            end
        end
        measure(t, SE, SC)
    end

    return average(SC)
end

function prun(S::StochPXSolver,
            p1::Dict{String,Vector{Any}},
            p2::Dict{String,Vector{Any}},
            MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
end

function average(SC::StochPXContext)
    nmesh = get_b("nmesh")
    ngrid = get_b("ngrid")
    ntry = get_x("ntry")

    Gout = zeros(C64, nmesh)
    Gr = zeros(F64, 2 * ngrid)
    c = 0.0
    for i = 1:ntry
        if SC.χ²[i] < 1e-6
            G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.mesh, SC.fmesh)
            @. Gout = Gout + G
            G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.grid, SC.fmesh)
            @. Gr = Gr + G
            c = c + 1.0
        end
    end
    @. Gout = Gout / c
    @. Gr = Gr / c

    return -imag.(Gout) / π, Gout, Gr
end

function last(SC::StochPXContext, Aout::Vector{F64}, Gout::Vector{C64}, Gr::Vector{F64})
    write_spectrum(SC.mesh, Aout)
    write_backward(SC.grid, Gr)
    write_complete(SC.mesh, Gout)
end

function warmup()
end

function sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    if rand(MC.rng) < 0.5
        try_move_p(t, MC, SE, SC)
    else
        try_move_a(t, MC, SE, SC)
    end
end

function measure(t::I64, SE::StochPXElement, SC::StochPXContext)
    @. SC.Pᵥ[t] = SE.P
    @. SC.Aᵥ[t] = SE.A
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

"""
    init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})

Randomize the configurations for future Monte Carlo sampling. It will
return a StochPXElement object.

See also: [`StochPXElement`](@ref).
"""
function init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})
    npole = get_x("npole")

    P = rand(rng, allow, npole)
    A = rand(rng, F64, npole)

    s = sum(A)
    @. A = A / s

    SE = StochPXElement(P, A)

    return SE
end

"""
    init_iodata(S::StochPXSolver, rd::RawData)

Preprocess the input data (`rd`), then allocate memory for the calculated
spectral functions.

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochPXSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹
end

function init_context(S::StochPXSolver)
    ntry = get_x("ntry")
    npole = get_x("npole")

    χ² = zeros(F64, ntry)

    Pᵥ = Vector{I64}[]
    Aᵥ = Vector{F64}[]
    for i = 1:ntry
        push!(Pᵥ, zeros(I64, npole))
        push!(Aᵥ, zeros(F64, npole))
    end

    return χ², Pᵥ, Aᵥ
end

function reset_mc(MC::StochPXMC)
    MC.Pacc = 0
    MC.Ptry = 0
    MC.Aacc = 0
    MC.Atry = 0
end

function reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
    npole = get_x("npole")

    P = rand(rng, allow, npole)
    A = rand(rng, F64, npole)

    s = sum(A)

    @. SE.P = P
    @. SE.A = A / s
end

function reset_iodata()
end

function reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
    ntry = get_x("ntry")
    @assert 1 ≤ t ≤ ntry

    Gᵧ = calc_green(SE.P, SE.A, SC.grid, SC.fmesh)
    χ² = calc_chi2(Gᵧ, SC.Gᵥ)
    @. SC.Gᵧ = Gᵧ
    SC.χ²[t] = χ²
end

"""
    calc_fmesh(S::StochPXSolver)

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to build the correlation function.

See also: [`LinearMesh`](@ref).
"""
function calc_fmesh(S::StochPXSolver)
    nfine = get_x("nfine")
    wmin = get_b("wmin")
    wmax = get_b("wmax")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
end

function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    grid::AbstractGrid,
                    fmesh::AbstractMesh)
    npole = get_x("npole")
    ngrid = length(grid)

    G = zeros(C64, ngrid)
    for i in eachindex(grid)
        z = 0.0
        for j = 1:npole
            z = z + A[j] / ( im * grid[i] - fmesh[P[j]] )
        end
        G[i] = z
    end

    return vcat(real(G), imag(G))
end

function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    mesh::AbstractMesh,
                    fmesh::AbstractMesh)
    npole = get_x("npole")
    nmesh = length(mesh)

    G = zeros(C64, nmesh)
    for i in eachindex(mesh)
        z = 0.0
        for j = 1:npole
            z = z + A[j] / ( mesh[i] + im * 1.0e-4 - fmesh[P[j]] )
        end
        G[i] = z
    end

    return G
end

function calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
    ΔG = Gₙ - Gᵥ
    return dot(ΔG, ΔG)
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

function try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    ngrid = get_b("ngrid")
    npole = get_x("npole")

    δG = zeros(C64, ngrid)

    s1 = 1
    s2 = 1
    while s1 == s2
        s1 = rand(MC.rng, 1:npole)
        s2 = rand(MC.rng, 1:npole)
    end

    P1 = SE.P[s1]
    P3 = P1
    while P3 == P1
        P3 = rand(MC.rng, SC.allow)
    end
    A1 = SE.A[s1]
    #
    P2 = SE.P[s2]
    P4 = P2
    while P4 == P2
        P4 = rand(MC.rng, SC.allow)
    end
    A2 = SE.A[s2]

    for i in eachindex(SC.grid)
        z = 0 + A1 / ( im * SC.grid[i] - SC.fmesh[P3] )
        z = z - A1 / ( im * SC.grid[i] - SC.fmesh[P1] )
        z = z + A2 / ( im * SC.grid[i] - SC.fmesh[P4] )
        z = z - A2 / ( im * SC.grid[i] - SC.fmesh[P2] )
        δG[i] = z
    end

    Gₙ = SC.Gᵧ + vcat(real(δG), imag(δG))
    χ² = calc_chi2(Gₙ, SC.Gᵥ)

    if χ² < SC.χ²[t]
        SE.P[s1] = P3
        SE.P[s2] = P4
        @. SC.Gᵧ = Gₙ
        SC.χ²[t] = χ²
    end
end

function try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    ngrid = get_b("ngrid")
    npole = get_x("npole")

    δG = zeros(C64, ngrid)

    s1 = 1
    s2 = 1
    while s1 == s2
        s1 = rand(MC.rng, 1:npole)
        s2 = rand(MC.rng, 1:npole)
    end

    P1 = SE.P[s1]
    P2 = SE.P[s2]
    A1 = SE.A[s1]
    A2 = SE.A[s2]
    A3 = 0
    A4 = 0
    while true
        δA = rand(MC.rng) * (A1 + A2) - A1
        A3 = A1 + δA
        A4 = A2 - δA

        if A3 > 0 && A4 > 0
            break
        end
    end

    for i in eachindex(SC.grid)
        z = 0 + (A3 - A1) / ( im * SC.grid[i] - SC.fmesh[P1] )
        z = z + (A4 - A2) / ( im * SC.grid[i] - SC.fmesh[P2] )
        δG[i] = z
    end

    Gₙ = SC.Gᵧ + vcat(real(δG), imag(δG))
    χ² = calc_chi2(Gₙ, SC.Gᵥ)

    if χ² < SC.χ²[t]
        SE.A[s1] = A3
        SE.A[s2] = A4
        @. SC.Gᵧ = Gₙ
        SC.χ²[t] = χ²
    end
end
