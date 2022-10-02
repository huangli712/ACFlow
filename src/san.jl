const P_SAC = Dict{String,Any}(
    "beta" => 8.0,
    "ntime" => 160,
    "nbins" => 1000,
    "nbootstrap" => 1000,
    "freq_interval" => 2.0e-5,
    "spec_interval" => 1.0e-2,
    "ommax" => 10.0,
    "ommin" => -10.0,
    "sac_bin_num" => 1,
    "sac_bin_size" => 100,
    "annealling_steps" => 1000,
    "ndelta" => 1000,
    "collecting_steps" => 1000,
    "stabilization_pace" => 10,
    "theta" => 1e+6,
    "annealing_rate" => 0.9
)

mutable struct SACContext
    Gr :: Vector{F64}
    G1 :: Vector{F64}
    G2 :: Vector{F64}
    χ2 :: F64
    χ2min :: F64
    Θ :: F64
    freq :: Vector{F64}
    spectrum :: Vector{F64}
end

struct SACGrid
    ommax :: F64
    ommin :: F64
    freq_interval :: F64
    spec_interval :: F64
    num_freq_index :: I64
    num_spec_index :: I64
end

mutable struct SACElement
    C :: Vector{I64}
    A :: F64
    W :: I64
end

struct SACAnnealing
    Conf  :: Vector{SACElement}
    Theta :: Vector{F64}
    chi2  :: Vector{F64}
end

mutable struct SACMonteCarlo
    rng :: AbstractRNG
    acc :: F64
    sample_acc  :: Vector{F64}
    sample_chi2 :: Vector{F64}
    bin_acc :: Vector{F64}
    bin_chi2 :: Vector{F64}
end

function FreqIndex2Freq(freq_index::I64, SG::SACGrid)
    @assert 1 ≤ freq_index ≤ SG.num_freq_index
    return SG.ommin + (freq_index - 1) * SG.freq_interval
end

function SpecIndex2Freq(spec_index::I64, SG::SACGrid)
    @assert 1 ≤ spec_index ≤ SG.num_spec_index
    return SG.ommin + (spec_index - 1) * SG.spec_interval
end

function calc_grid()
    ommax = P_SAC["ommax"]
    ommin = P_SAC["ommin"]
    freq_interval = P_SAC["freq_interval"]
    spec_interval = P_SAC["spec_interval"]
    num_freq_index = ceil(I64, (ommax - ommin) / freq_interval)
    num_spec_index = ceil(I64, (ommax - ommin) / spec_interval)

    return SACGrid(ommax, ommin, freq_interval, spec_interval, num_freq_index, num_spec_index)
end

function read_data()
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

function init_kernel(tmesh, SG::SACGrid, Mrot::AbstractMatrix)
    beta = P_SAC["beta"]

    ntau = length(tmesh)
    nfreq = SG.num_freq_index
    kernel = zeros(F64, ntau, nfreq)

    for f = 1:nfreq
        ω = FreqIndex2Freq(f, SG)
        de = 1.0 + exp(-beta * ω)
        kernel[:,f] = exp.(-ω * tmesh) / de
    end

    kernel = Mrot * kernel

    return kernel
end

function init_mc()
    sbin = P_SAC["sac_bin_size"]
    nbin = P_SAC["sac_bin_num"]

    seed = rand(1:1000000); seed = 840443
    rng = MersenneTwister(seed)
    acc = 0.0
    sample_acc = zeros(F64, sbin)
    sample_chi2 = zeros(F64, sbin)
    bin_acc = zeros(F64, nbin)
    bin_chi2 = zeros(F64, nbin)
    MC = SACMonteCarlo(rng, acc, sample_acc, sample_chi2, bin_acc, bin_chi2)

    return MC
end

function init_spectrum(rng, scale_factor::F64, SG::SACGrid, Gdata, tau)
    ndelta = P_SAC["ndelta"]

    position = zeros(I64, ndelta)
    rand!(rng, position, 1:SG.num_freq_index)

    amplitude = 1.0 / (scale_factor * ndelta)
    average_freq = abs(log(1.0/Gdata[end]) / tau[end])
    window_width = ceil(I64, 0.1 * average_freq / SG.freq_interval)

    return SACElement(position, amplitude, window_width)
end

function compute_corr_from_spec(kernel::AbstractMatrix, SE::SACElement, SC::SACContext)
    ndelta = P_SAC["ndelta"]
    tmp_kernel = kernel[:, SE.C]
    amplitude = fill(SE.A, ndelta)
    SC.G1 = tmp_kernel * amplitude
end

function compute_goodness(G::Vector{F64,}, Gr::Vector{F64}, Sigma::Vector{F64})
    χ = sum(((G .- Gr) .* Sigma) .^ 2.0)
    return χ
end

function perform_annealing(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, covar)
    anneal_length = P_SAC["annealling_steps"]

    Conf = SACElement[]
    Theta = F64[]
    Chi2 = F64[]

    for i = 1:anneal_length
        update_fixed_theta(MC, SE, SC, SG, kernel, covar)

        SC.χ2 = mean(MC.bin_chi2)

        push!(Conf, deepcopy(SE))
        push!(Theta, SC.Θ)
        push!(Chi2, SC.χ2)

        @show i, SC.χ2, SC.χ2min, SC.χ2 - SC.χ2min
        if SC.χ2 - SC.χ2min < 1e-3
            break
        end

        SC.Θ = SC.Θ * P_SAC["annealing_rate"]
    end

    return SACAnnealing(Conf, Theta, Chi2)
end

function decide_sampling_theta(anneal::SACAnnealing, SC::SACContext, kernel::AbstractMatrix, covar)
    num_anneal = length(anneal.chi2)

    c = num_anneal
    while c ≥ 1
        if anneal.chi2[c] > SC.χ2min + 2.0 * sqrt(SC.χ2min)
            break
        end
        c = c - 1
    end
    @assert 1 ≤ c ≤ num_anneal

    SE = deepcopy(anneal.Conf[c])
    SC.Θ = anneal.Theta[c]
    compute_corr_from_spec(kernel, SE, SC)
    SC.χ2 = compute_goodness(SC.G1, SC.Gr, covar)
    @show SC.Θ, SC.χ2

    return SE
end

function update_fixed_theta(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, covar)
    nbin = P_SAC["sac_bin_num"]
    sbin = P_SAC["sac_bin_size"]
    ntau = length(covar)

    for n = 1:nbin
        for s = 1:sbin

            if (s - 1) % P_SAC["stabilization_pace"] == 1
                SC.χ2 = compute_goodness(SC.G1, SC.Gr, covar)
            end

            update_deltas_1step_single(MC, SE, SC, SG, kernel, covar)

            MC.sample_chi2[s] = SC.χ2
            MC.sample_acc[s] = MC.acc
        end

        MC.bin_chi2[n] = sum(MC.sample_chi2) / sbin
        MC.bin_acc[n] = sum(MC.sample_acc) / sbin

        @show n, SC.Θ, SC.χ2min / ntau, MC.bin_chi2[n] / ntau,  MC.bin_chi2[n] - SC.χ2min, MC.bin_acc[n], SE.W * SG.freq_interval

        if MC.bin_acc[n] > 0.5
            r = SE.W * 1.5
            if ceil(I64, r) < SG.num_freq_index
                SE.W = ceil(I64, r)
            else
                SE.W = SG.num_freq_index
            end
        end

        if MC.bin_acc[n] < 0.4
            SE.W = ceil(I64, SE.W / 1.5)
        end
    end
end

function update_deltas_1step_single(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, covar)
    ndelta = P_SAC["ndelta"]
    accept_count = 0.0

    for i = 1:ndelta
        select_delta = rand(MC.rng, 1:ndelta)
        location_current = SE.C[select_delta]

        if 1 < SE.W < SG.num_freq_index
            move_width = rand(MC.rng, 1:SE.W)

            if rand(MC.rng) > 0.5
                location_updated = location_current + move_width
            else
                location_updated = location_current - move_width
            end

            if location_updated < 1 
                location_updated = location_updated + SG.num_freq_index
            end

            if location_updated > SG.num_freq_index
                location_updated = location_updated - SG.num_freq_index
            end

        elseif SE.W == SG.num_freq_index
            location_updated = rand(MC.rng, 1:SG.num_freq_index)
        else
            error("BIG PROBLEM")
        end

        SC.G2 = SC.G1 + SE.A .* (kernel[:,location_updated] .- kernel[:,location_current])

        chi2_updated = compute_goodness(SC.G2, SC.Gr, covar)

        p = exp( (SC.χ2 - chi2_updated) / (2.0 * SC.Θ) )

        if rand(MC.rng) < min(p, 1.0)
            SE.C[select_delta] = location_updated
            SC.G1 = deepcopy(SC.G2)
            SC.χ2 = chi2_updated
            if SC.χ2 < SC.χ2min
                SC.χ2min = SC.χ2
            end

            accept_count = accept_count + 1.0
        end
    end

    MC.acc = accept_count / ndelta
end

function sample_and_collect(scale_factor::F64, MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::AbstractMatrix, covar)
    ndelta = P_SAC["ndelta"]
    update_fixed_theta(MC, SE, SC, SG, kernel, covar)

    for n = 1:SG.num_spec_index
        SC.freq[n] = SpecIndex2Freq(n, SG)
    end

    collecting_steps = P_SAC["collecting_steps"]
    for i = 1:collecting_steps
        if (i - 1) % P_SAC["stabilization_pace"] == 1
            SC.χ2 = compute_goodness(SC.G1, SC.Gr, covar)
            @show i, SC.χ2
        end

        update_deltas_1step_single(MC, SE, SC, SG, kernel, covar)

        for j = 1:ndelta
            d_pos = SE.C[j]
            s_pos = Grid2Spec(d_pos, SG)
            #@show d_pos, s_pos
            SC.spectrum[s_pos] = SC.spectrum[s_pos] + SE.A
        end
    end

    factor = scale_factor / (collecting_steps * SG.spec_interval)
    SC.spectrum = SC.spectrum * factor

    open("Aout.data", "w") do fout
        for i = 1:SG.num_spec_index
            println(fout, SC.freq[i], " ", SC.spectrum[i])
        end
    end
end

function Grid2Spec(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_freq_index
    return ceil(I64, grid_index * SG.freq_interval / SG.spec_interval)
end

function sac_run()
    tmesh, gbin = read_data()
    factor, gtau = compute_corr_means(gbin)
    gerr, bootstrape = compute_corr_errs(gbin, gtau)
    tmesh, gerr, gtau, bootstrape = discard_poor_quality_data(tmesh, gerr, gtau, bootstrape)
    gtau, gerr, bootstrape = scale_data(factor, gtau, gerr, bootstrape)
    vals, vecs, cov_mat = compute_cov_matrix(gtau, bootstrape)
    covar = calc_covar(vals)

    grid = calc_grid()
    kernel = init_kernel(tmesh, grid, vecs)
    mc = init_mc()
    SE = init_spectrum(mc.rng, factor, grid, gtau, tmesh)

    ntau = length(tmesh)
    Gr = vecs * gtau
    G1 = zeros(F64, ntau)
    G2 = zeros(F64, ntau)
    χ2 = 0.0
    χ2min = 0.0
    Θ = P_SAC["theta"]
    freq = zeros(F64, grid.num_spec_index)
    spectrum = zeros(F64, grid.num_spec_index)
    SC = SACContext(Gr, G1, G2, χ2, χ2min, Θ, freq, spectrum)
    compute_corr_from_spec(kernel, SE, SC)
    χ = compute_goodness(SC.G1, SC.Gr, covar)
    SC.χ2 = χ
    SC.χ2min = χ
    anneal = perform_annealing(mc, SE, SC, grid, kernel, covar)
    SE = decide_sampling_theta(anneal, SC, kernel, covar)
    sample_and_collect(factor, mc, SE, SC, grid, kernel, covar)
end
