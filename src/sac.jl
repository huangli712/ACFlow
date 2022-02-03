module StochAC

using Random
using LinearAlgebra

import ..ACFlow: I64, F64, C64
import ..ACFlow: AbstractMesh, UniformMesh, AbstractGrid
import ..ACFlow: RawData
import ..ACFlow: make_data, make_grid, make_mesh, make_model, make_kernel
import ..ACFlow: get_c, get_a

#=
const P_Stoch = Dict{String,Any}(
    "ngrid" => 1000,
    "nmesh" => 801,
    "nfine" => 10001,
    "ngamm" => 1024,
    "nalph" => 6,
    "nwarm" => 4000,
    "nstep" => 4000000,
    "ndump" => 40000,
    "alpha" => 1.00,
    "ratio" => 2.00,
    "beta"  => 5.00,
)
=#

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
    alist  :: Vector{F64}
    hamil  :: Vector{F64}
    HC     :: Array{F64,2}
    tmesh :: AbstractGrid
    wmesh :: AbstractMesh
    G_tau :: Vector{F64}
    G_dev :: Vector{F64}
end

mutable struct StochMC
    rng :: AbstractRNG
    move_acc :: Vector{I64}
    move_try :: Vector{I64}
    swap_acc :: Vector{I64}
    swap_try :: Vector{I64}
end

function solve(rd::RawData)
    G = make_data(rd)
    Gtau = abs.(G.value) ./ G.var
    Gdev = G.var
    tmesh = make_grid(rd)

    MC, SE, SC = stoch_init(tmesh, Gtau, Gdev)
    stoch_run(MC, SE, SC)
end

function stoch_grid()
    nfine = get_a("nfine")
    wmin = get_c("wmin")
    wmax = get_c("wmax")

    fmesh = UniformMesh(nfine, wmin, wmax)

    model = fill(1.0/nfine, nfine)
    xmesh = cumsum(model)

    return fmesh, xmesh
end

function stoch_delta(xmesh::Vector{F64}, phi::Vector{F64})
    nmesh = get_c("nmesh")
    nfine = get_a("nfine")
    eta1 = 0.005
    eta2 = 0.005 ^ 2.0

    delta = zeros(F64, nmesh, nfine)
    
    for i = 1:nfine
        s = phi .- xmesh[i]
        delta[:,i] = eta1 ./ (s .* s .+ eta2)
    end

    return delta
end

function stoch_init(tmesh::AbstractGrid, G_tau::Vector{F64}, G_dev::Vector{F64})
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")
    alpha = get_a("alpha")
    ratio = get_a("ratio")
    nfine = get_a("nfine")
    ngamm = get_a("ngamm")
    ngrid = get_c("ngrid")

    seed = rand(1:100000000); seed = 4277216
    @show seed
    rng = MersenneTwister(seed)
    move_acc = zeros(F64, nalph)
    move_try = zeros(F64, nalph)
    swap_acc = zeros(F64, nalph)
    swap_try = zeros(F64, nalph)
    MC = StochMC(rng, move_acc, move_try, swap_acc, swap_try)

    r_γ = rand(rng, F64, (ngamm, nalph))
    a_γ = rand(rng, 1:nfine, (ngamm, nalph))
    for j = 1:nalph
        s = sum(r_γ[:,j])
        r_γ[:,j] = r_γ[:,j] ./ s
    end
    SE = StochElement(a_γ, r_γ)

    wmesh = make_mesh()
    fmesh, xmesh = stoch_grid()

    alist = collect(alpha * (ratio ^ (x - 1)) for x in 1:nalph)
    model = make_model(wmesh)
    phi = cumsum(model .* wmesh.weight)
    delta = stoch_delta(xmesh, phi)

    kernel = make_kernel(fmesh, tmesh)

    image = zeros(F64, nmesh, nalph)
    hamil = zeros(F64, nalph)
    HC = zeros(F64, ngrid, nalph)
    δt = tmesh[2] - tmesh[1]
    for i = 1:nalph
        HC[:,i] = stoch_hamil0(a_γ[:,i], r_γ[:,i], kernel, G_tau, G_dev)
        hamil[i] = dot(HC[:,i], HC[:,i]) * δt
    end
    SC = StochContext(kernel, delta, image, phi, model, alist, hamil, HC, tmesh, wmesh, G_tau, G_dev)

    return MC, SE, SC
end

function stoch_hamil0(a_γ::Vector{I64},
                      r_γ::Vector{F64},
                      kernel::Array{F64,2},
                      G_tau::Vector{F64},
                      G_dev::Vector{F64})
    ngrid = get_c("ngrid")
    ngamm = get_a("ngamm")

    hc = zeros(F64, ngrid)
    for i = 1:ngrid
        hc[i] = dot(r_γ, kernel[i,a_γ])
    end
    @. hc = hc / G_dev - G_tau

    return hc
end

function stoch_dump(step::F64, MC::StochMC, SC::StochContext)
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
            println(fout, "# $i :", SC.alist[i])
            for j = 1:nmesh
                println(fout, SC.wmesh[j], " ", image_t[j,i])
            end
        end
    end

    open("stoch.data.sum", "w") do fout
        for j = 1:nmesh
            println(fout, SC.wmesh[j], " ", sum(image_t[j,:]) / nalph)
        end
    end
end

function stoch_run(MC::StochMC, SE::StochElement, SC::StochContext)
    nstep = get_a("nstep")
    ndump = get_a("ndump")

    stoch_warmming(MC, SE, SC)

    step = 0.0 
    for iter = 1:nstep
        stoch_sampling(MC, SE, SC)
            
        if iter % 100 == 0
            step = step + 1.0
            stoch_recording(SE, SC)
        end

        if iter % ndump == 0
            println("iter: ", iter / ndump)
            stoch_dump(step, MC, SC)
        end
    end
end

function stoch_recording(SE::StochElement, SC::StochContext)
    nalph = get_a("nalph")
    nmesh = get_c("nmesh")
    ngamm = get_a("ngamm")

    for ia = 1:nalph
        dr = view(SE.r_γ, :, ia)
        da = view(SE.a_γ, :, ia)
        Aw = SC.delta[:,da] * dr
        SC.image[:,ia] = SC.image[:,ia] .+ Aw
    end
end

function try_mov1(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    ngamm = get_a("ngamm")

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
    dhc = dhh * (K1 - K2) ./ SC.G_dev

    δt = SC.tmesh[2] - SC.tmesh[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 ||  exp(-SC.alist[i] * dhh) > rand(MC.rng)
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

function try_mov2(i::I64, MC::StochMC, SE::StochElement, SC::StochContext)
    ngamm = get_a("ngamm")
    nfine = get_a("nfine")

    hc = view(SC.HC, :, i)

    l1 = 1
    l2 = 1
    while l1 == l2
        l1 = rand(MC.rng, 1:ngamm)
        l2 = rand(MC.rng, 1:ngamm)
    end

    r1 = SE.r_γ[l1,i]
    r2 = SE.r_γ[l2,i]

    i1 = rand(MC.rng, 1:nfine)
    i2 = rand(MC.rng, 1:nfine)
    i3 = SE.a_γ[l1,i]
    i4 = SE.a_γ[l2,i]

    K1 = view(SC.kernel, :, i1)
    K2 = view(SC.kernel, :, i2)
    K3 = view(SC.kernel, :, i3)
    K4 = view(SC.kernel, :, i4)
    dhc = ( r1 * (K1 - K3) + r2 * (K2 - K4) ) ./ SC.G_dev

    δt = SC.tmesh[2] - SC.tmesh[1]
    dhh = dot(dhc, 2.0 * hc + dhc) * δt

    pass = false
    if dhh ≤ 0.0 ||  exp(-SC.alist[i] * dhh) > rand(MC.rng)
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

    da = SC.alist[i] - SC.alist[j]
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

function stoch_sampling(MC::StochMC, SE::StochElement, SC::StochContext)
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

function stoch_warmming(MC::StochMC, SE::StochElement, SC::StochContext)
    nwarm = get_a("nwarm")

    for i = 1:nwarm
        println("warm: $i")
        stoch_sampling(MC, SE, SC)
    end

    fill!(MC.move_acc, 0.0)
    fill!(MC.move_try, 0.0)

    fill!(MC.swap_acc, 0.0)
    fill!(MC.swap_try, 0.0)
end
end