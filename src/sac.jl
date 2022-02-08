#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/08
#

"""
    StochMC
"""
mutable struct StochMC
    rng :: AbstractRNG
    move_acc :: Vector{I64}
    move_try :: Vector{I64}
    swap_acc :: Vector{I64}
    swap_try :: Vector{I64}
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
    σ²     :: Vector{F64}
    grid   :: AbstractGrid
    mesh   :: AbstractMesh
    model  :: Vector{F64}
    kernel :: Array{F64,2}
    image  :: Array{F64,2}
    Δ      :: Array{F64,2}
    hτ     :: Array{F64,2}
    Hα     :: Vector{F64}
    αₗ     :: Vector{F64}
end

"""
    solve
"""
function solve(S::StochACSolver, rd::RawData)
    println("[ StochAC ]")
    MC, SE, SC = init(S, rd)
    run(S, MC, SE, SC)
end

"""
    init
"""
function init(S::StochACSolver, rd::RawData)
    MC = init_mc()
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(MC.rng)
    println("Randomize Monte Carlo configurations")

    Gᵥ, σ², image = init_iodata(rd)
    println("Postprocess input data: ", length(σ²), " points")

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

    hτ, Hα = calc_hamil(SE, grid, kernel, Gᵥ, σ²)
    println("Precompute hamiltonian")

    αₗ = calc_alpha()
    println("Precompute α parameters")

    SC = StochContext(Gᵥ, σ², grid, mesh, model, kernel, image, Δ, hτ, Hα, αₗ)

    return MC, SE, SC
end

"""
    run
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

"""
    warmup
"""
function warmup(MC::StochMC, SE::StochElement, SC::StochContext)
    nwarm = get_a("nwarm")

    for i = 1:nwarm
        #println("warm: $i")
        sample(MC, SE, SC)
    end

    fill!(MC.move_acc, 0.0)
    fill!(MC.move_try, 0.0)

    fill!(MC.swap_acc, 0.0)
    fill!(MC.swap_try, 0.0)
end

"""
    sample
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
    measure
"""
function measure(SE::StochElement, SC::StochContext)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")
    ngamm = get_a("ngamm")

    for ia = 1:nalph
        dr = view(SE.Γᵣ, :, ia)
        da = view(SE.Γₐ, :, ia)
        Aw = SC.Δ[:,da] * dr
        SC.image[:,ia] = SC.image[:,ia] .+ Aw
    end
end

"""
    init_mc()

Try to create a StochMC struct.

See also: [`StochAC`](@ref).
"""
function init_mc()
    nalph = get_a("nalph")

    seed = rand(1:100000000)

    rng = MersenneTwister(seed)
    move_acc = zeros(F64, nalph)
    move_try = zeros(F64, nalph)
    swap_acc = zeros(F64, nalph)
    swap_try = zeros(F64, nalph)

    MC = StochMC(rng, move_acc, move_try, swap_acc, swap_try)

    return MC
end

"""
    init_element()
"""
function init_element(rng::AbstractRNG)
    nalph = get_a("nalph")
    nfine = get_a("nfine")
    ngamm = get_a("ngamm")

    Γᵣ = rand(rng, F64, (ngamm, nalph))
    Γₐ = rand(rng, 1:nfine, (ngamm, nalph))

    for j = 1:nalph
        s = sum(Γᵣ[:,j])
        Γᵣ[:,j] = Γᵣ[:,j] ./ s
    end

    SE = StochElement(Γₐ, Γᵣ)

    return SE
end

"""
    init_iodata(rd::RawData)
"""
function init_iodata(rd::RawData)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")

    G = make_data(rd)
    Gᵥ = abs.(G.value)
    σ² = 1.0 ./ G.covar

    image = zeros(F64, nmesh, nalph)

    return Gᵥ, σ², image
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
    calc_phi(mesh::AbstractMesh, model::Vector{F64})

Try to calculate ϕ(ω) function.

See also: [`AbstractMesh`](@ref).
"""
function calc_phi(mesh::AbstractMesh, model::Vector{F64})
    ϕ = cumsum(model .* mesh.weight)
    return ϕ
end

"""
    calc_delta()
"""
function calc_delta(xmesh::Vector{F64}, ϕ::Vector{F64})
    nmesh = get_c("nmesh")
    nfine = get_a("nfine")
    η₁ = 0.005
    η₂ = 0.005 ^ 2.0

    Δ = zeros(F64, nmesh, nfine)

    for i = 1:nfine
        s = ϕ .- xmesh[i]
        Δ[:,i] = η₁ ./ (s .* s .+ η₂)
    end

    return Δ
end

"""
    calc_hamil(...)

Initialize h(τ) and H using Eq.(35) and Eq.(36).
"""
function calc_hamil(SE::StochElement, grid::AbstractGrid, kernel, Gᵥ, σ²)
    nalph = get_a("nalph")
    ngrid = get_c("ngrid")

    hτ = zeros(F64, ngrid, nalph)
    Hα = zeros(F64, nalph)

    δt = grid[2] - grid[1]
    for i = 1:nalph
        hτ[:,i] = calc_htau(SE.Γₐ[:,i], SE.Γᵣ[:,i], kernel, Gᵥ, σ²)
        Hα[i] = dot(hτ[:,i], hτ[:,i]) * δt
    end

    return hτ, Hα
end

"""
    calc_htau(...)

Try to calculate h(τ) via Eq.(36).
"""
function calc_htau(Γₐ::Vector{I64}, Γᵣ::Vector{F64},
                   kernel::Matrix{F64},
                   Gᵥ::Vector{F64}, σ²::Vector{F64})
    ngrid = get_c("ngrid")

    hτ = zeros(F64, ngrid)
    for i = 1:ngrid
        hτ[i] = dot(Γᵣ, kernel[i,Γₐ])
    end
    @. hτ = (hτ - Gᵥ) * σ²

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
    try_mov1()
"""
function try_mov1(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    ngamm = get_a("ngamm")

    hc = view(SC.hτ, :, i)

    l1 = 1
    l2 = 1
    while l1 == l2
        l1 = rand(MC.rng, 1:ngamm)
        l2 = rand(MC.rng, 1:ngamm)
    end

    r1 = 0.0
    r2 = 0.0
    dhh = 0.0
    r3 = SE.Γᵣ[l1,i]
    r4 = SE.Γᵣ[l2,i]
    while true
        dhh = rand(MC.rng) * (r3 + r4) - r3
        r1 = r3 + dhh
        r2 = r4 - dhh
        if r1 > 0 && r2 > 0
            break
        end
    end

    i1 = SE.Γₐ[l1,i]
    i2 = SE.Γₐ[l2,i]

    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    dhc = dhh * (K1 - K2) .* SC.σ²

    δt = SC.grid[2] - SC.grid[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 || exp(-SC.αₗ[i] * dhh) > rand(MC.rng)
        pass = true
    end

    if pass
        SE.Γᵣ[l1,i] = r1
        SE.Γᵣ[l2,i] = r2
        @. hc = hc + dhc
        SC.Hα[i] = dot(hc, hc) * δt
    end

    MC.move_try[i] = MC.move_try[i] + 1.0
    if pass
        MC.move_acc[i] = MC.move_acc[i] + 1.0
    end
end

"""
    try_mov2()
"""
function try_mov2(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    ngamm = get_a("ngamm")
    nfine = get_a("nfine")

    hc = view(SC.hτ, :, i)

    l1 = 1
    l2 = 1
    while l1 == l2
        l1 = rand(MC.rng, 1:ngamm)
        l2 = rand(MC.rng, 1:ngamm)
    end

    r1 = SE.Γᵣ[l1,i]
    r2 = SE.Γᵣ[l2,i]

    i1 = rand(MC.rng, 1:nfine)
    i2 = rand(MC.rng, 1:nfine)
    i3 = SE.Γₐ[l1,i]
    i4 = SE.Γₐ[l2,i]

    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    K3 = view(SC.kernel, :, i3)
    K4 = view(SC.kernel, :, i4)
    dhc = ( r1 * (K1 - K3) + r2 * (K2 - K4) ) .* SC.σ²

    δt = SC.grid[2] - SC.grid[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 ||  exp(-SC.αₗ[i] * dhh) > rand(MC.rng)
        pass = true
    end

    if pass
        SE.Γₐ[l1,i] = i1
        SE.Γₐ[l2,i] = i2
        @. hc = hc + dhc
        SC.Hα[i] = dot(hc, hc) * δt
    end

    MC.move_try[i] = MC.move_try[i] + 1.0
    if pass
        MC.move_acc[i] = MC.move_acc[i] + 1.0
    end
end

"""
    try_swap()
"""
function try_swap(scheme::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    nalph = get_a("nalph")

    if scheme == 1
        i = rand(MC.rng, 1:nalph)
        if rand(MC.rng) > 0.5
            j = i + 1
        else
            j = i - 1
        end

        i == 1 && (j = i + 1)
        i == nalph && (j = i - 1)
    else
        while true
            i = rand(MC.rng, 1:nalph)
            j = rand(MC.rng, 1:nalph)
            i != j && break
        end
    end

    da = SC.αₗ[i] - SC.αₗ[j]
    dh = SC.Hα[i] - SC.Hα[j]

    pass = ( exp(da * dh) > rand(MC.rng) )

    if pass
        SE.Γₐ[:,i], SE.Γₐ[:,j] = SE.Γₐ[:,j], SE.Γₐ[:,i]
        SE.Γᵣ[:,i], SE.Γᵣ[:,j] = SE.Γᵣ[:,j], SE.Γᵣ[:,i]

        SC.hτ[:,i], SC.hτ[:,j] = SC.hτ[:,j], SC.hτ[:,i]
        SC.Hα[i], SC.Hα[j] = SC.Hα[j], SC.Hα[i]
    end

    MC.swap_try[i] = MC.swap_try[i] + 1.0
    if pass
        MC.swap_acc[i] = MC.swap_acc[i] + 1.0
    end
end

function dump(step::F64, MC::StochMC, SC::StochContext)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")

    println("move statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.move_acc[i] / MC.move_try[i])
    end
    println("swap statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.swap_acc[i] / MC.swap_try[i])
    end

    image_t = zeros(F64, nmesh, nalph)
    for i = 1:nalph
        for j = 1:nmesh
            image_t[j,i] = SC.image[j,i] * SC.model[j] / π / step
        end
    end

    open("stoch.data", "w") do fout
        for i = 1:nalph
            println(fout, "# $i :", SC.αₗ[i])
            for j = 1:nmesh
                println(fout, SC.mesh[j], " ", image_t[j,i])
            end
        end
    end

    open("stoch.data.sum", "w") do fout
        for j = 1:nmesh
            println(fout, SC.mesh[j], " ", sum(image_t[j,:]) / nalph)
        end
    end
end
