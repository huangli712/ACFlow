#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/20
#

#=
### *Customized Structs* : *StochAC Solver*
=#

"""
    StochMC
"""
mutable struct StochMC
    rng :: AbstractRNG
    Macc :: Vector{I64}
    Mtry :: Vector{I64}
    Sacc :: Vector{I64}
    Stry :: Vector{I64}
end

"""
    StochElement
"""
mutable struct StochElement
    Γₐ :: Array{I64,2}
    Γᵣ :: Array{F64,2}
end

"""
    StochContext
"""
mutable struct StochContext
    Gᵥ     :: Vector{F64}
    σ¹     :: Vector{F64}
    grid   :: AbstractGrid
    mesh   :: AbstractMesh
    model  :: Vector{F64}
    kernel :: Array{F64,2}
    Aout   :: Array{F64,2}
    Δ      :: Array{F64,2}
    hτ     :: Array{F64,2}
    Hα     :: Vector{F64}
    αₗ     :: Vector{F64}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochACSolver, rd::RawData)
"""
function solve(S::StochACSolver, rd::RawData)
    println("[ StochAC ]")
    MC, SE, SC = init(S, rd)
    run(S, MC, SE, SC)
end

"""
    init(S::StochACSolver, rd::RawData)
"""
function init(S::StochACSolver, rd::RawData)
    MC = init_mc()
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(MC.rng)
    println("Randomize Monte Carlo configurations")

    Gᵥ, σ¹, Aout = init_iodata(rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    model = make_model(mesh)
    println("Build default model: ", get_c("mtype"))

    fmesh = calc_fmesh()
    kernel = make_kernel(fmesh, grid)
    println("Build default kernel: ", get_c("ktype"))

    xmesh = calc_xmesh()
    ϕ = calc_phi(mesh, model)
    Δ = calc_delta(xmesh, ϕ)
    println("Precompute δ functions")

    hτ, Hα = calc_hamil(SE.Γₐ, SE.Γᵣ, grid, kernel, Gᵥ, σ¹)
    println("Precompute hamiltonian")

    αₗ = calc_alpha()
    println("Precompute α parameters")

    SC = StochContext(Gᵥ, σ¹, grid, mesh, model, kernel, Aout, Δ, hτ, Hα, αₗ)

    return MC, SE, SC
end

"""
    run(S::StochACSolver, MC::StochMC, SE::StochElement, SC::StochContext)
"""
function run(S::StochACSolver, MC::StochMC, SE::StochElement, SC::StochContext)
    nstep = get_a("nstep")
    ndump = get_a("ndump")

    warmup(MC, SE, SC)

    step = 0.0
    for iter = 1:nstep
        sample(MC, SE, SC)

        if iter % 100 == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % ndump == 0
            println("iter: ", iter / ndump)
            dump(step, MC, SC)
        end
    end
end

function postprocess()
end

#=
### *Core Algorithms*
=#

"""
    warmup(MC::StochMC, SE::StochElement, SC::StochContext)
"""
function warmup(MC::StochMC, SE::StochElement, SC::StochContext)
    nwarm = get_a("nwarm")

    for i = 1:nwarm
        sample(MC, SE, SC)
    end

    fill!(MC.Macc, 0.0)
    fill!(MC.Mtry, 0.0)

    fill!(MC.Sacc, 0.0)
    fill!(MC.Stry, 0.0)
end

"""
    sample(MC::StochMC, SE::StochElement, SC::StochContext)
"""
function sample(MC::StochMC, SE::StochElement, SC::StochContext)
    nalph = get_a("nalph")

    if rand(MC.rng) < 0.9
        if rand(MC.rng) > 0.5
            for i = 1:nalph
                try_mov1(i, MC, SE, SC)
            end
        else
            for i = 1:nalph
                try_mov2(i, MC, SE, SC)
            end
        end
    else
        if nalph > 1
            if rand(MC.rng) > 0.1
                try_swap(1, MC, SE, SC)
            else
                try_swap(2, MC, SE, SC)
            end
        end
    end
end

"""
    measure(SE::StochElement, SC::StochContext)
"""
function measure(SE::StochElement, SC::StochContext)
    nalph = get_a("nalph")

    for ia = 1:nalph
        dr = view(SE.Γᵣ, :, ia)
        da = view(SE.Γₐ, :, ia)
        Aw = SC.Δ[:,da] * dr
        SC.Aout[:,ia] = SC.Aout[:,ia] .+ Aw
    end
end

#=
### *Service Functions*
=#

"""
    init_mc()

Try to create a StochMC struct.

See also: [`StochAC`](@ref).
"""
function init_mc()
    nalph = get_a("nalph")

    seed = rand(1:100000000); seed = 39061530
    @show seed

    rng = MersenneTwister(seed)
    Macc = zeros(F64, nalph)
    Mtry = zeros(F64, nalph)
    Sacc = zeros(F64, nalph)
    Stry = zeros(F64, nalph)

    MC = StochMC(rng, Macc, Mtry, Sacc, Stry)

    return MC
end

"""
    init_element(rng::AbstractRNG)

Randomize the configurations for future monte carlo sampling. It will
return a StochElement object.

See also: [`StochElement`](@ref).
"""
function init_element(rng::AbstractRNG)
    nalph = get_a("nalph")
    nfine = get_a("nfine")
    ngamm = get_a("ngamm")

    Γᵣ = rand(rng, F64, (ngamm, nalph))
    Γₐ = rand(rng, 1:nfine, (ngamm, nalph))

    for j = 1:nalph
        Γⱼ = view(Γᵣ, :, j)
        s = sum(Γⱼ)
        @. Γⱼ = Γⱼ / s
    end

    SE = StochElement(Γₐ, Γᵣ)

    return SE
end

"""
    init_iodata(rd::RawData)

Preprocess the input data (`rd`), then allocate memory for the α-resolved
spectral functions.

See also: [`RawData`](@ref).
"""
function init_iodata(rd::RawData)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")

    Aout = zeros(F64, nmesh, nalph)

    G = make_data(rd)
    Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹, Aout
end

"""
    calc_fmesh()

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to build the kernel function.

See also: [`LinearMesh`](@ref).
"""
function calc_fmesh()
    nfine = get_a("nfine")
    wmin = get_c("wmin")
    wmax = get_c("wmax")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
end

"""
    calc_xmesh()

Try to calculate very fine (dense) linear mesh in [0, 1], which is used
internally to build the δ functions.

See also: [`calc_delta`](@ref).
"""
function calc_xmesh()
    nfine = get_a("nfine")

    _mesh = fill(1.0/nfine, nfine)
    xmesh = cumsum(_mesh)

    return xmesh
end

"""
    calc_phi(am::AbstractMesh, model::Vector{F64})

Try to calculate ϕ(ω) function. `am` is the mesh for calculated spectrum,
and `model` means the default model function.

See also: [`AbstractMesh`](@ref), [`calc_delta`](@ref).
"""
function calc_phi(am::AbstractMesh, model::Vector{F64})
    ϕ = cumsum(model .* am.weight)
    return ϕ
end

"""
    calc_delta(xmesh::Vector{F64}, ϕ::Vector{F64})

Precompute the δ functions. `xmesh` is a very dense linear mesh in [0,1]
and `ϕ` is the ϕ function.

See also: [`calc_xmesh`](@ref), [`calc_phi`](@ref).å
"""
function calc_delta(xmesh::Vector{F64}, ϕ::Vector{F64})
    nmesh = length(ϕ)
    nfine = length(xmesh)

    η₁ = 0.001
    η₂ = 0.001 ^ 2.0

    Δ = zeros(F64, nmesh, nfine)
    s = similar(ϕ)
    for i = 1:nfine
        @. s = (ϕ - xmesh[i]) ^ 2.0 + η₂
        @. Δ[:,i] = η₁ / s
    end

    return Δ
end

"""
    calc_hamil(Γₐ, Γᵣ, grid, kernel, Gᵥ, σ¹)

Initialize h(τ) and H(α) using Eq.(35) and Eq.(36), respectively. `Γₐ`
and `Γᵣ` represent n(x), `grid` is the grid for the input data, `kernel`
means the kernel function, `Gᵥ` is the correlator, and `σ¹` is equal to
1.0 / σ.
"""
function calc_hamil(Γₐ::Array{I64,2}, Γᵣ::Array{F64,2},
                    grid::AbstractGrid,
                    kernel::Matrix{F64},
                    Gᵥ::Vector{F64}, σ¹::Vector{F64})
    nalph = get_a("nalph")
    ngrid = length(Gᵥ)

    hτ = zeros(F64, ngrid, nalph)
    Hα = zeros(F64, nalph)

    δt = grid[2] - grid[1]
    for i = 1:nalph
        hτ[:,i] = calc_htau(Γₐ[:,i], Γᵣ[:,i], kernel, Gᵥ, σ¹)
        Hα[i] = dot(hτ[:,i], hτ[:,i]) * δt
    end

    return hτ, Hα
end

"""
    calc_htau(Γₐ, Γᵣ, kernel, Gᵥ, σ¹)

Try to calculate α-dependent h(τ) via Eq.(36). `Γₐ` and `Γᵣ` represent
n(x), `kernel` means the kernel function, `Gᵥ` is the correlator, and
`σ¹` is equal to 1.0 / σ.

See also: [`calc_hamil`](@ref).
"""
function calc_htau(Γₐ::Vector{I64}, Γᵣ::Vector{F64}, 
                   kernel::Matrix{F64},
                   Gᵥ::Vector{F64}, σ¹::Vector{F64})
    hτ = similar(Gᵥ)

    for i in eachindex(Gᵥ)
        hτ[i] = dot(Γᵣ, view(kernel, i, Γₐ))
    end

    @. hτ = (hτ - Gᵥ) * σ¹

    return hτ
end

"""
    calc_alpha()

Generate a list for the α parameters
"""
function calc_alpha()
    nalph = get_a("nalph")
    alpha = get_a("alpha")
    ratio = get_a("ratio")

    αₗ = collect(alpha * (ratio ^ (x - 1)) for x in 1:nalph)

    return αₗ
end

"""
    try_mov1(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)

Select two δ functions and then change their weights. Here `i` means the
index for α parameters.

See also: [`try_mov2`](@ref).
"""
function try_mov1(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    # Get current number of δ functions
    ngamm = get_a("ngamm")

    # Choose two δ functions, they are labelled as γ1 and γ2, respectively.
    γ1 = 1
    γ2 = 1
    while γ1 == γ2
        γ1 = rand(MC.rng, 1:ngamm)
        γ2 = rand(MC.rng, 1:ngamm)
    end

    # Extract weights for the two δ functions (r3 and r4), then try to
    # calculate new weights for them (r1 and r2).
    r1 = 0.0
    r2 = 0.0
    r3 = SE.Γᵣ[γ1,i]
    r4 = SE.Γᵣ[γ2,i]
    δr = 0.0
    while true
        δr = rand(MC.rng) * (r3 + r4) - r3
        r1 = r3 + δr
        r2 = r4 - δr
        if r1 > 0 && r2 > 0
            break
        end
    end

    # Try to calculate the change of Hc using Eq.~(42).
    hc = view(SC.hτ, :, i)
    K1 = view(SC.kernel, :, SE.Γₐ[γ1,i])
    K2 = view(SC.kernel, :, SE.Γₐ[γ2,i])
    #
    δt = SC.grid[2] - SC.grid[1]
    δhc = δr * (K1 - K2) .* SC.σ¹
    δH = dot(δhc, 2.0 * hc + δhc) * δt

    # Apply Metropolis algorithm
    MC.Mtry[i] = MC.Mtry[i] + 1.0
    if δH ≤ 0.0 || exp(-SC.αₗ[i] * δH) > rand(MC.rng)
        # Update monte carlo configurations
        SE.Γᵣ[γ1,i] = r1
        SE.Γᵣ[γ2,i] = r2

        # Update h(τ)
        @. hc = hc + δhc

        # Update Hc
        SC.Hα[i] = SC.Hα[i] + δH

        # Update monte carlo counter
        MC.Macc[i] = MC.Macc[i] + 1.0
    end
end

"""
    try_mov2(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)

Select two δ functions and then change their positions. Here `i` means the
index for α parameters.

See also: [`try_mov1`](@ref).
"""
function try_mov2(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    # Get current number of δ functions
    ngamm = get_a("ngamm")

    # Get total number of δ functions
    nfine = get_a("nfine")

    # Choose two δ functions, they are labelled as γ1 and γ2, respectively.
    γ1 = 1
    γ2 = 1
    while γ1 == γ2
        γ1 = rand(MC.rng, 1:ngamm)
        γ2 = rand(MC.rng, 1:ngamm)
    end

    # Extract weights for the two δ functions (r1 and r2)
    r1 = SE.Γᵣ[γ1,i]
    r2 = SE.Γᵣ[γ2,i]

    # Choose new positions for the two δ functions (i1 and i2).
    # Note that their old positions are SE.Γₐ[γ1,i] and SE.Γₐ[γ2,i].
    i1 = rand(MC.rng, 1:nfine)
    i2 = rand(MC.rng, 1:nfine)

    # Try to calculate the change of Hc using Eq.~(42).
    hc = view(SC.hτ, :, i)
    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    K3 = view(SC.kernel, :, SE.Γₐ[γ1,i])
    K4 = view(SC.kernel, :, SE.Γₐ[γ2,i])
    #
    δt = SC.grid[2] - SC.grid[1]
    δhc = ( r1 * (K1 - K3) + r2 * (K2 - K4) ) .* SC.σ¹
    δH = dot(δhc, 2.0 * hc + δhc) * δt

    # Apply Metropolis algorithm
    MC.Mtry[i] = MC.Mtry[i] + 1.0
    if δH ≤ 0.0 || exp(-SC.αₗ[i] * δH) > rand(MC.rng)
        # Update monte carlo configurations
        SE.Γₐ[γ1,i] = i1
        SE.Γₐ[γ2,i] = i2

        # Update h(τ)
        @. hc = hc + δhc

        # Update Hc
        SC.Hα[i] = SC.Hα[i] + δH

        # Update monte carlo counter
        MC.Macc[i] = MC.Macc[i] + 1.0
    end
end

"""
    try_swap(MC::StochMC, SE::StochElement, SC::StochContext)

Try to exchange field configurations between two adjacent layers.
"""
function try_swap(MC::StochMC, SE::StochElement, SC::StochContext)
    # Get number of α parameters
    nalph = get_a("nalph")

    # Select two adjacent layers (two adjacent α parameters)
    i = rand(MC.rng, 1:nalph)
    j = rand(MC.rng) > 0.5 ? i + 1 : i - 1
    i == 1 && (j = i + 1)
    i == nalph && (j = i - 1)
 
    # Calculate change of Hc
    δα = SC.αₗ[i] - SC.αₗ[j]
    δH = SC.Hα[i] - SC.Hα[j]

    # Apply Metropolis algorithm
    MC.Stry[i] = MC.Stry[i] + 1.0
    MC.Stry[j] = MC.Stry[j] + 1.0
    if exp(δα * δH) > rand(MC.rng)
        # Update monte carlo configurations
        SE.Γₐ[:,i], SE.Γₐ[:,j] = SE.Γₐ[:,j], SE.Γₐ[:,i]
        SE.Γᵣ[:,i], SE.Γᵣ[:,j] = SE.Γᵣ[:,j], SE.Γᵣ[:,i]

        # Update h(τ) and Hc
        SC.hτ[:,i], SC.hτ[:,j] = SC.hτ[:,j], SC.hτ[:,i]
        SC.Hα[i], SC.Hα[j] = SC.Hα[j], SC.Hα[i]

        # Update monte carlo counter
        MC.Sacc[i] = MC.Sacc[i] + 1.0
        MC.Sacc[j] = MC.Sacc[j] + 1.0
    end
end

function dump(step::F64, MC::StochMC, SC::StochContext)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")

    println("move statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.Macc[i] / MC.Mtry[i])
    end
    println("swap statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.Sacc[i] / MC.Stry[i])
    end

    Aw = zeros(F64, nmesh, nalph)
    for i = 1:nalph
        for j = 1:nmesh
            Aw[j,i] = SC.Aout[j,i] * SC.model[j] / π / step
        end
    end
    Asum = [sum(Aw[i,:]) / nalph for i = 1:nmesh]

    write_hamil(SC.αₗ, SC.Hα)
    write_spectrum(SC.mesh, SC.αₗ, Aw)
    write_spectrum(SC.mesh, Asum)
end
