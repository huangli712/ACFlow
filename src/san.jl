const P_SAC = Dict{String,Any}(
    "ntime" => 160,
    "nbins" => 1000,
    "nbootstrap" => 1000,
)

mutable struct StochSKElement
    P :: Vector{I64}
    A :: F64
    W :: I64
end

mutable struct StochSKContext
    Gᵥ     :: Vector{F64}
    Gᵧ     :: Vector{F64}
    σ¹     :: Vector{F64}
    mesh   :: AbstractMesh
    kernel :: Array{F64,2}
    Aout   :: Vector{F64}
    χ²     :: F64
    χ²min  :: F64
    χ²vec  :: Vector{F64}
    Θ      :: F64
    Θvec   :: Vector{F64}
    𝒞ᵧ     :: Vector{StochSKElement}
end

function read_gtau()
    nbins = P_SAC["nbins"]
    ntime = P_SAC["ntime"]

    tgrid = zeros(F64, ntime)
    open("tgrids.in", "r") do fin
        readline(fin)
        for i = 1:ntime
            arr = line_to_array(fin)
            tgrid[i] = parse(F64, arr[2])
        end
    end

    gbin = zeros(F64, nbins, ntime)
    open("corr.in", "r") do fin
        readline(fin)
        for i = 1:nbins
            for j = 1:ntime
                arr = line_to_array(fin)
                gbin[i,j] = parse(F64, arr[3])
            end
        end
    end

    return tgrid, gbin
end

function compute_corr_means(gbin)
    A = vec(mean(gbin, dims = 1))
    factor = A[1]
    return factor, A
end

function compute_corr_errs(gbin, gtau)
    nbins = P_SAC["nbins"]
    ntime = P_SAC["ntime"]
    nbootstrap = P_SAC["nbootstrap"]

    gerr = zeros(F64, ntime)
    bootstrap_samples = zeros(F64, nbootstrap, ntime)

    rng = MersenneTwister(rand(1:10000) + 1981)
    for i = 1:nbootstrap
        v = zeros(F64, ntime)
        for _ = 1:nbins
            k = rand(rng, 1:nbins)
            v = v + gbin[k,:]
        end
        bootstrap_samples[i,:] = v[:]
    end
    bootstrap_samples = bootstrap_samples ./ nbins

    for i = 1:ntime
        for j = 1:nbootstrap
            gerr[i] = gerr[i] + (bootstrap_samples[j,i] - gtau[i]) ^ 2.0
        end
        gerr[i] = sqrt(gerr[i] / nbootstrap)
    end

    return gerr, bootstrap_samples
end

function discard_poor_quality_data(tmesh, gerr, gtau, bootstrap_samples)
    ntime = P_SAC["ntime"]
    good_tgrids = I64[]
    for i = 2:ntime
        if abs(gerr[i] / gtau[i]) < 0.1
            push!(good_tgrids, i)
        end
    end

    tmesh = tmesh[good_tgrids]
    gtau = gtau[good_tgrids]
    gerr = gerr[good_tgrids]
    bootstrap_samples = bootstrap_samples[:, good_tgrids]

    return tmesh, gerr, gtau, bootstrap_samples
end

function scale_data(factor, gtau, gerr, bootstrape)
    gtau = gtau ./ factor
    gerr = gerr ./ factor
    bootstrape = bootstrape ./ factor

    return gtau, gerr, bootstrape
end

function calc_covar(vals)
    nbootstrap = P_SAC["nbootstrap"]
    covar = zeros(F64, length(vals))
    for i in eachindex(vals)
        covar[i] = sqrt(nbootstrap) / sqrt(vals[i])
    end
    return covar
end

function compute_cov_matrix(gtau, bootstrap_samples)
    ncov = length(gtau)
    cov_mat = zeros(F64, ncov, ncov)

    for i = 1:ncov
        for j = 1:ncov
            cov_mat[i,j] = sum((bootstrap_samples[:,i] .- gtau[i]) .* (bootstrap_samples[:,j] .- gtau[j]))
        end
    end

    F = eigen(cov_mat)
    vals, vecs = F

    return vals, vecs, cov_mat
end

function san_run()
    tmesh, gbin = read_gtau()
    factor, gtau = compute_corr_means(gbin)
    gerr, bootstrape = compute_corr_errs(gbin, gtau)
    tmesh, gerr, gtau, bootstrape = discard_poor_quality_data(tmesh, gerr, gtau, bootstrape)
    gtau, gerr, bootstrape = scale_data(factor, gtau, gerr, bootstrape)
    vals, vecs, cov_mat = compute_cov_matrix(gtau, bootstrape)
    σ¹ = calc_covar(vals)

    fmesh = LinearMesh(get_k("nfine"), get_b("wmin"), get_b("wmax"))
    kernel = init_kernel(tmesh, fmesh, vecs)
    mc = init_mc()
    SE = init_delta(mc.rng, factor, fmesh, gtau, tmesh)

    ntau = length(tmesh)
    Gᵥ = vecs * gtau
    Gᵧ = zeros(F64, ntau)
    χ² = 0.0
    χ²min = 0.0
    χ²vec = zeros(F64, get_k("nwarm"))
    Θ = get_k("theta")
    Θvec = zeros(F64, get_k("nwarm"))
    mesh = LinearMesh(get_b("nmesh"), get_b("wmin"), get_b("wmax"))
    Aout = zeros(F64, get_b("nmesh"))
    𝒞ᵧ = StochSKElement[]
    SC = StochSKContext(Gᵥ, Gᵧ, σ¹, mesh, kernel, Aout, χ², χ²min, χ²vec, Θ, Θvec, 𝒞ᵧ)
    compute_corr_from_spec(SE, SC)
    χ = compute_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
    SC.χ² = χ
    SC.χ²min = χ
    warmup(mc, SE, SC, fmesh)
    SE = decide_sampling_theta(SC)
    measure(factor, mc, SE, SC, fmesh)
end

function init_kernel(tmesh, fmesh::AbstractMesh, Mrot::AbstractMatrix)
    beta = get_b("beta")

    ntau = length(tmesh)
    nfreq = length(fmesh)
    kernel = zeros(F64, ntau, nfreq)

    for f = 1:nfreq
        ω = fmesh[f]
        de = 1.0 + exp(-beta * ω)
        kernel[:,f] = exp.(-ω * tmesh) / de
    end

    kernel = Mrot * kernel

    return kernel
end

function init_mc()
    seed = rand(1:1000000); seed = 840443
    rng = MersenneTwister(seed)
    acc = 0.0
    MC = StochSKMC(rng, acc)

    return MC
end

function init_delta(rng, scale_factor::F64, fmesh::AbstractMesh, Gdata, tau)
    ngamm = get_k("ngamm")

    position = zeros(I64, ngamm)
    rand!(rng, position, 1:length(fmesh))

    amplitude = 1.0 / (scale_factor * ngamm)
    average_freq = abs(log(1.0/Gdata[end]) / tau[end])
    window_width = ceil(I64, 0.1 * average_freq / (fmesh[2] - fmesh[1]))

    return StochSKElement(position, amplitude, window_width)
end

function warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext, fmesh::AbstractMesh)
    anneal_length = get_k("nwarm")

    for i = 1:anneal_length
        SC.χ² = update_fixed_theta(MC, SE, SC, fmesh)

        push!(SC.𝒞ᵧ, deepcopy(SE))
        SC.χ²vec[i] = SC.χ²
        SC.Θvec[i] = SC.Θ

        @show i, SC.χ², SC.χ²min, SC.χ² - SC.χ²min
        if SC.χ² - SC.χ²min < 1e-3
            break
        end

        SC.Θ = SC.Θ * get_k("ratio")
    end
end

function decide_sampling_theta(SC::StochSKContext)
    num_anneal = length(SC.𝒞ᵧ)
    @assert num_anneal ≤ get_k("nwarm")

    c = num_anneal
    while c ≥ 1
        if SC.χ²vec[c] > SC.χ²min + 2.0 * sqrt(SC.χ²min)
            break
        end
        c = c - 1
    end
    @assert 1 ≤ c ≤ num_anneal

    SE = deepcopy(SC.𝒞ᵧ[c])
    SC.Θ = SC.Θvec[c]
    compute_corr_from_spec(SE, SC)
    SC.χ² = compute_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
    @show SC.Θ, SC.χ²

    return SE
end

function measure(scale_factor::F64, MC::StochSKMC, SE::StochSKElement, SC::StochSKContext, fmesh::AbstractMesh)
    ngamm = get_k("ngamm")
    update_fixed_theta(MC, SE, SC, fmesh)

    nstep = get_k("nstep")
    for i = 1:nstep
        if (i - 1) % get_k("retry") == 1
            SC.χ² = compute_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
            @show i, SC.χ²
        end

        update_deltas_1step_single(MC, SE, SC, fmesh)

        for j = 1:ngamm
            d_pos = SE.P[j]
            s_pos = ceil(I64, d_pos / length(fmesh) * get_b("nmesh"))
            SC.Aout[s_pos] = SC.Aout[s_pos] + SE.A
        end
    end

    factor = scale_factor / (nstep * (SC.mesh[2] - SC.mesh[1]))
    SC.Aout = SC.Aout * factor

    open("Aout.data", "w") do fout
        for i in eachindex(SC.mesh)
            println(fout, SC.mesh[i], " ", SC.Aout[i])
        end
    end
end

function compute_corr_from_spec(SE::StochSKElement, SC::StochSKContext)
    ngamm = get_k("ngamm")
    tmp_kernel = SC.kernel[:, SE.P]
    amplitude = fill(SE.A, ngamm)
    SC.Gᵧ = tmp_kernel * amplitude
end

function compute_goodness(G::Vector{F64,}, Gᵥ::Vector{F64}, Sigma::Vector{F64})
    χ = sum(((G .- Gᵥ) .* Sigma) .^ 2.0)
    return χ
end

function update_fixed_theta(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext, fmesh::AbstractMesh)
    nbin = 1
    sbin = 100
    ntau = length(SC.σ¹)

    sample_chi2 = zeros(F64, sbin)
    bin_chi2 = zeros(F64, nbin)
    sample_acc = zeros(F64, sbin)
    bin_acc = zeros(F64, nbin)

    for n = 1:nbin
        for s = 1:sbin

            if (s - 1) % get_k("retry") == 1
                SC.χ² = compute_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
            end

            update_deltas_1step_single(MC, SE, SC, fmesh)

            sample_chi2[s] = SC.χ²
            sample_acc[s] = MC.acc
        end

        bin_chi2[n] = sum(sample_chi2) / sbin
        bin_acc[n] = sum(sample_acc) / sbin

        @show n, SC.Θ, SC.χ²min / ntau, bin_chi2[n] / ntau,  bin_chi2[n] - SC.χ²min, bin_acc[n], SE.W * (fmesh[2] - fmesh[1])

        if bin_acc[n] > 0.5
            r = SE.W * 1.5
            if ceil(I64, r) < length(fmesh)
                SE.W = ceil(I64, r)
            else
                SE.W = length(fmesh)
            end
        end

        if bin_acc[n] < 0.4
            SE.W = ceil(I64, SE.W / 1.5)
        end
    end

    return mean(bin_chi2)
end

function update_deltas_1step_single(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext, fmesh::AbstractMesh)
    ngamm = get_k("ngamm")
    accept_count = 0.0

    for i = 1:ngamm
        select_delta = rand(MC.rng, 1:ngamm)
        location_current = SE.P[select_delta]

        if 1 < SE.W < length(fmesh)
            move_width = rand(MC.rng, 1:SE.W)

            if rand(MC.rng) > 0.5
                location_updated = location_current + move_width
            else
                location_updated = location_current - move_width
            end

            if location_updated < 1 
                location_updated = location_updated + length(fmesh)
            end

            if location_updated > length(fmesh)
                location_updated = location_updated - length(fmesh)
            end

        elseif SE.W == length(fmesh)
            location_updated = rand(MC.rng, 1:length(fmesh))
        else
            error("BIG PROBLEM")
        end

        Gₙ = SC.Gᵧ + SE.A .* (SC.kernel[:,location_updated] .- SC.kernel[:,location_current])

        chi2_updated = compute_goodness(Gₙ, SC.Gᵥ, SC.σ¹)

        p = exp( (SC.χ² - chi2_updated) / (2.0 * SC.Θ) )

        if rand(MC.rng) < min(p, 1.0)
            SE.P[select_delta] = location_updated
            SC.Gᵧ = deepcopy(Gₙ)
            SC.χ² = chi2_updated
            if SC.χ² < SC.χ²min
                SC.χ²min = SC.χ²
            end

            accept_count = accept_count + 1.0
        end
    end

    MC.acc = accept_count / ngamm
end
