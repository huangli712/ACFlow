module StochAC

using Random
using LinearAlgebra

const I64 = Int64
const F64 = Float64
const C64 = ComplexF64

export read_data
export stoch_init
export stoch_run

const P_Stoch = Dict{String,Any}(
    "ntime" => 1000,
    "nwmax" => 801,
    "ngrid" => 10001,
    "ngamm" => 1024,
    "nalph" => 6,
    "nwarm" => 4000,
    "nstep" => 4000000,
    "ndump" => 40000,
    "ainit" => 1.00,
    "ratio" => 2.00,
    "beta"  => 5.00,
    "eta1"  => 0.005,
    "eta2"  => 0.005^2,
    "sigma" => 0.0,
    "wstep" => 0.02
)

mutable struct GreenData
    G_qmc :: Vector{F64}
    G_tau :: Vector{F64}
    G_dev :: Vector{F64}
end

mutable struct StochGrid
    tmesh :: Vector{F64}
    wmesh :: Vector{F64}
    xgrid :: Vector{F64}
    wgrid :: Vector{F64}
end

mutable struct StochElement
    a_γ :: Array{I64,2}
    r_γ :: Array{F64,2}
end

mutable struct StochContext
    kernel :: Array{F64,2}
    delta  :: Array{F64,2}
    image  :: Array{F64,2}
    phi    :: Vector{F64}
    model  :: Vector{F64}
    alpha  :: Vector{F64}
    hamil  :: Vector{F64}
    HC     :: Array{F64,2}
end

mutable struct StochMC
    rng :: AbstractRNG
    move_acc :: Vector{I64}
    move_try :: Vector{I64}
    swap_acc :: Vector{I64}
    swap_try :: Vector{I64}
end

function read_data()
    ntime = P_Stoch["ntime"]
    G_qmc = zeros(F64, ntime)
    G_tau = zeros(F64, ntime)
    G_dev = zeros(F64, ntime)
    tmesh = zeros(F64, ntime)

#=    
    open("solver.green.dat", "r") do fin
        for i = 1:ntime
            arr = split(readline(fin))
            tmesh[i] = parse(F64, arr[3])
            G_qmc[i] = parse(F64, arr[4])
            G_dev[i] = parse(F64, arr[5])
            if abs(G_dev[i]) < 1e-6
                G_dev[i] = 1e-6
            end
            G_tau[i] = abs(G_qmc[i]) / G_dev[i]
        end
    end
=#

    open("green.data", "r") do fin
        for i = 1:ntime
            arr = split(readline(fin))
            tmesh[i] = parse(F64, arr[1])
            G_qmc[i] = parse(F64, arr[2])
            G_dev[i] = parse(F64, arr[3])
            if abs(G_dev[i]) < 1e-6
                G_dev[i] = 1e-6
            end
            G_dev[i] = G_dev[i] ^ 2.0
            G_tau[i] = abs(G_qmc[i]) / G_dev[i]
        end
    end

    return GreenData(G_qmc, G_tau, G_dev), tmesh
end

function stoch_norm!(weight::F64, fun::AbstractVector{F64})
    norm = sum(fun) * weight
    @. fun = fun / norm
end

function stoch_grid()
    ngrid = P_Stoch["ngrid"]
    nwmax = P_Stoch["nwmax"]
    wstep = P_Stoch["wstep"]

    ommin = -(nwmax - 1) / 2 * wstep
    ommax = +(nwmax - 1) / 2 * wstep
    wgrid = collect(LinRange(ommin, ommax, ngrid))

    model = fill(1.0, ngrid)
    stoch_norm!(1.0, model)
    xgrid = cumsum(model)

    return wgrid, xgrid
end

function stoch_delta(xgrid::Vector{F64}, phi::Vector{F64})
    nwmax = P_Stoch["nwmax"]
    ngrid = P_Stoch["ngrid"]
    eta1 = P_Stoch["eta1"]
    eta2 = P_Stoch["eta2"]

    delta = zeros(F64, nwmax, ngrid)
    
    for i = 1:ngrid
        s = phi .- xgrid[i]
        delta[:,i] = eta1 ./ (s .* s .+ eta2)
    end

    return delta
end

function stoch_kernel(tmesh::Vector{F64}, wgrid::Vector{F64})
    ntime = P_Stoch["ntime"]
    ngrid = P_Stoch["ngrid"]
    β = P_Stoch["beta"]

    kernel = zeros(F64, ntime, ngrid)

    for j = 1:ngrid
        ω = wgrid[j]
        if ω ≥ 0.0
            denom = 1.0 + exp(-β*ω)
            kernel[:,j] = exp.(-tmesh * ω) / denom
        else
            denom = 1.0 + exp(+β*ω)
            kernel[:,j] = exp.( ( β .- tmesh ) * ω) / denom
        end
    end
   
    return kernel
end

function stoch_init(tmesh::Vector{F64}, G::GreenData)
    nalph = P_Stoch["nalph"]
    nwmax = P_Stoch["nwmax"]
    wstep = P_Stoch["wstep"]
    ainit = P_Stoch["ainit"]
    ratio = P_Stoch["ratio"]
    ngrid = P_Stoch["ngrid"]
    ngamm = P_Stoch["ngamm"]
    ntime = P_Stoch["ntime"]

    seed = rand(1:100000000)#; seed = 19087549
    @show seed
    rng = MersenneTwister(seed)
    move_acc = zeros(F64, nalph)
    move_try = zeros(F64, nalph)
    swap_acc = zeros(F64, nalph)
    swap_try = zeros(F64, nalph)
    MC = StochMC(rng, move_acc, move_try, swap_acc, swap_try)

    r_γ = rand(rng, F64, (ngamm, nalph))
    a_γ = rand(rng, 1:ngrid, (ngamm, nalph))
    for j = 1:nalph
        stoch_norm!(1.0, view(r_γ, :, j))
    end
    SE = StochElement(a_γ, r_γ)

    ommin = -(nwmax - 1) / 2 * wstep
    ommax = +(nwmax - 1) / 2 * wstep
    wmesh = collect(LinRange(ommin, ommax, nwmax))
    wgrid, xgrid = stoch_grid()
    SG = StochGrid(tmesh, wmesh, xgrid, wgrid)

    model = fill(1.0, nwmax)
    stoch_norm!(wstep, model)
    alpha = collect(ainit * (ratio ^ (x - 1)) for x in 1:nalph)
    phi = cumsum(model) * wstep
    delta = stoch_delta(xgrid, phi)
    kernel = stoch_kernel(tmesh, wgrid)
    image = zeros(F64, nwmax, nalph)
    hamil = zeros(F64, nalph)
    HC = zeros(F64, ntime, nalph)
    δt = tmesh[2] - tmesh[1]
    for i = 1:nalph
        HC[:,i] = stoch_hamil0(a_γ[:,i], r_γ[:,i], kernel, G.G_tau, G.G_dev)
        hamil[i] = dot(HC[:,i], HC[:,i]) * δt
    end
    SC = StochContext(kernel, delta, image, phi, model, alpha, hamil, HC)

    return MC, SE, SG, SC
end

function stoch_hamil0(a_γ::Vector{I64},
                      r_γ::Vector{F64},
                      kernel::Array{F64,2},
                      G_tau::Vector{F64},
                      G_dev::Vector{F64})
    ntime = P_Stoch["ntime"]
    ngamm = P_Stoch["ngamm"]

    hc = zeros(F64, ntime)
    for i = 1:ntime
        hc[i] = dot(r_γ, kernel[i,a_γ])
    end
    @. hc = hc / G_dev - G_tau

    return hc
end

function stoch_dump(step::F64, MC::StochMC, SC::StochContext, SG::StochGrid)
    nalph = P_Stoch["nalph"]
    nwmax = P_Stoch["nwmax"]

    println("move statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.move_acc[i] / MC.move_try[i])
    end
    println("swap statistics:")
    for i = 1:nalph
        println("    alpha $i: ", MC.swap_acc[i] / MC.swap_try[i])
    end

    image_t = zeros(F64, nwmax, nalph)
    for i = 1:nalph
        for j = 1:nwmax
            image_t[j,i] = SC.image[j,i] * SC.model[j] / π / step
        end
    end

    open("stoch.data", "w") do fout
        for i = 1:nalph
            println(fout, "# $i :", SC.alpha[i])
            for j = 1:nwmax
                println(fout, SG.wmesh[j], " ", image_t[j,i])
            end
        end
    end

    open("stoch.data.sum", "w") do fout
        for j = 1:nwmax
            println(fout, SG.wmesh[j], " ", sum(image_t[j,:]) / nalph)
        end
    end
end

function stoch_run(MC::StochMC, SE::StochElement, SC::StochContext, SG::StochGrid, G::GreenData)
    nstep = P_Stoch["nstep"]
    ndump = P_Stoch["ndump"]

    stoch_warmming(MC, SE, SC, SG, G)

    step = 0.0 
    for iter = 1:nstep
        stoch_sampling(MC, SE, SC, SG, G)
            
        if iter % 100 == 0
            step = step + 1.0
            stoch_recording(SE, SC)
        end

        if iter % ndump == 0
            println("iter: ", iter / ndump)
            stoch_dump(step, MC, SC, SG)
        end
    end
end

function stoch_recording(SE::StochElement, SC::StochContext)
    nalph = P_Stoch["nalph"]
    nwmax = P_Stoch["nwmax"]
    ngamm = P_Stoch["ngamm"]

    for ia = 1:nalph
        dr = view(SE.r_γ, :, ia)
        da = view(SE.a_γ, :, ia)
        Aw = SC.delta[:,da] * dr
        SC.image[:,ia] = SC.image[:,ia] .+ Aw
    end
end

function try_mov1(i::I64, MC::StochMC, SE::StochElement, SC::StochContext, SG::StochGrid, G::GreenData)
    ngamm = P_Stoch["ngamm"]

    hc = view(SC.HC, :, i)

    l1 = 1
    l2 = 1
    while l1 == l2
        l1 = rand(MC.rng, 1:ngamm)
        l2 = rand(MC.rng, 1:ngamm)
    end

    r1 = 0.0
    r2 = 0.0
    dhh = 0.0
    r3 = SE.r_γ[l1,i]
    r4 = SE.r_γ[l2,i]
    while true
        dhh = rand(MC.rng) * (r3 + r4) - r3
        r1 = r3 + dhh
        r2 = r4 - dhh
        if r1 > 0 && r2 > 0
            break
        end
    end

    i1 = SE.a_γ[l1,i]
    i2 = SE.a_γ[l2,i]

    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    dhc = dhh * (K1 - K2) ./ G.G_dev

    δt = SG.tmesh[2] - SG.tmesh[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 ||  exp(-SC.alpha[i] * dhh) > rand(MC.rng)
        pass = true
    end

    if pass
        SE.r_γ[l1,i] = r1
        SE.r_γ[l2,i] = r2
        @. hc = hc + dhc
        SC.hamil[i] = dot(hc, hc) * δt
    end

    MC.move_try[i] = MC.move_try[i] + 1.0
    if pass
        MC.move_acc[i] = MC.move_acc[i] + 1.0
    end
end

function try_mov2(i::I64, MC::StochMC, SE::StochElement, SC::StochContext, SG::StochGrid, G::GreenData)
    ngamm = P_Stoch["ngamm"]
    ngrid = P_Stoch["ngrid"]

    hc = view(SC.HC, :, i)

    l1 = 1
    l2 = 1
    while l1 == l2
        l1 = rand(MC.rng, 1:ngamm)
        l2 = rand(MC.rng, 1:ngamm)
    end

    r1 = SE.r_γ[l1,i]
    r2 = SE.r_γ[l2,i]

    i1 = rand(MC.rng, 1:ngrid)
    i2 = rand(MC.rng, 1:ngrid)
    i3 = SE.a_γ[l1,i]
    i4 = SE.a_γ[l2,i]

    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    K3 = view(SC.kernel, :, i3)
    K4 = view(SC.kernel, :, i4)
    dhc = ( r1 * (K1 - K3) + r2 * (K2 - K4) ) ./ G.G_dev

    δt = SG.tmesh[2] - SG.tmesh[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 ||  exp(-SC.alpha[i] * dhh) > rand(MC.rng)
        pass = true
    end

    if pass
        SE.a_γ[l1,i] = i1
        SE.a_γ[l2,i] = i2
        @. hc = hc + dhc
        SC.hamil[i] = dot(hc, hc) * δt
    end

    MC.move_try[i] = MC.move_try[i] + 1.0
    if pass
        MC.move_acc[i] = MC.move_acc[i] + 1.0
    end
end

function try_swap(scheme::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    nalph = P_Stoch["nalph"]

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

    da = SC.alpha[i] - SC.alpha[j]
    dh = SC.hamil[i] - SC.hamil[j]

    pass = ( exp(da * dh) > rand(MC.rng) )

    if pass        
        SE.a_γ[:,i], SE.a_γ[:,j] = SE.a_γ[:,j], SE.a_γ[:,i]
        SE.r_γ[:,i], SE.r_γ[:,j] = SE.r_γ[:,j], SE.r_γ[:,i]
        
        SC.HC[:,i], SC.HC[:,j] = SC.HC[:,j], SC.HC[:,i]
        SC.hamil[i], SC.hamil[j] = SC.hamil[j], SC.hamil[i]
    end

    MC.swap_try[i] = MC.swap_try[i] + 1.0
    if pass
        MC.swap_acc[i] = MC.swap_acc[i] + 1.0
    end
end

function stoch_sampling(MC::StochMC, SE::StochElement, SC::StochContext, SG::StochGrid, G::GreenData)
    nalph = P_Stoch["nalph"]

    if rand(MC.rng) < 0.9
        if rand(MC.rng) > 0.5
            for i = 1:nalph
                try_mov1(i, MC, SE, SC, SG, G)
            end
        else
            for i = 1:nalph
                try_mov2(i, MC, SE, SC, SG, G)
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

function stoch_warmming(MC::StochMC, SE::StochElement, SC::StochContext, SG::StochGrid, G::GreenData)
    nwarm = P_Stoch["nwarm"]

    for i = 1:nwarm
        println("warm: $i")
        stoch_sampling(MC, SE, SC, SG, G)
    end

    fill!(MC.move_acc, 0.0)
    fill!(MC.move_try, 0.0)

    fill!(MC.swap_acc, 0.0)
    fill!(MC.swap_try, 0.0)
end
end