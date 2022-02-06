#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/02/06
#

"""
    MaxEntContext
"""
mutable struct MaxEntContext
    Gᵥ     :: Vector{F64}
    σ²     :: Vector{F64}
    mesh   :: AbstractMesh
    model  :: Vector{F64}
    kernel :: Array{F64,2}
    hess   :: Array{F64,2}
    Vₛ     :: Array{F64,2}
    W₂     :: Array{F64,2}
    W₃     :: Array{F64,3}
    Bₘ     :: Vector{F64}
end

"""
    solve
"""
function solve(S::MaxEntSolver, rd::RawData)
    println("[ MaxEnt ]")
    mec = init(S, rd)
    run(S, mec)
end

"""
    init
"""
function init(S::MaxEntSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value
    σ² = 1.0 ./ G.covar
    println("Postprocess input data: ", length(σ²), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    model = make_model(mesh)
    println("Build default model: ", get_c("mtype"))

    kernel = make_kernel(mesh, grid)
    println("Build default kernel: ", get_c("ktype"))

    Vₛ, W₂, W₃, Bₘ, hess = precompute(Gᵥ, σ², mesh, model, kernel)
    println("Precompute key coefficients")

    return MaxEntContext(Gᵥ, σ², mesh, model, kernel, hess, Vₛ, W₂, W₃, Bₘ)
end

"""
    run
"""
function run(S::MaxEntSolver, mec::MaxEntContext)
    method = get_m("method")

    @cswitch method begin
        @case "historic"
            darr, sol = historic(mec)
            break

        @case "classic"
            darr, sol = classic(mec)
            break

        @case "bryan"
            darr, sol = bryan(mec)
            break

        @case "chi2kink"
            darr, sol = chi2kink(mec)
            break
    end

    postprocess(S, mec.mesh, darr, sol)
end

"""
    postprocess
"""
function postprocess(S::MaxEntSolver, am::AbstractMesh, darr, sol)
    write_spectrum(am, sol[:A_opt])
end

"""
    historic
"""
function historic(mec::MaxEntContext)
    println("Apply historic algorithm to determine optimized α")

    use_bayes = false
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)

    ustart = zeros(F64, n_svd)
    optarr = []

    conv = 0.0
    while conv < 1.0
        sol = optimizer(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / ratio
        conv = length(mec.σ²) / sol[:chi2]
    end

    ustart = optarr[end-1][:u_opt]
    alpha = optarr[end][:alpha]

    function root_fun(_alpha, _u)
        res = optimizer(mec, _alpha, _u, use_bayes)
        @. _u = res[:u_opt]
        return length(mec.σ²) / res[:chi2] - 1.0
    end
    alpha_opt = secant(root_fun, alpha, ustart)
    println("opt alpha:", alpha_opt)

    sol = optimizer(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol
end

"""
    classic
"""
function classic(mec::MaxEntContext)
    println("Using classic algorithm to solve the maximum entropy problem")

    use_bayes = true
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    conv = 0.0
    while conv < 1.0
        sol = optimizer(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / ratio
        @. ustart = sol[:u_opt]
        conv = sol[:conv]
    end

    convarr = [x[:conv] for x in optarr]
    alpharr = [x[:alpha] for x in optarr]
    exp_opt = log10(alpharr[end] / alpharr[end-1])
    exp_opt = exp_opt / log10(convarr[end] / convarr[end-1])
    exp_opt = log10(alpharr[end-1]) - log10(convarr[end-1]) * exp_opt
    alpha = 10.0 ^ exp_opt
    ustart = optarr[end-1][:u_opt]

    function root_fun(_alpha, _u)
        res = optimizer(mec, _alpha, _u, use_bayes)
        @. _u = res[:u_opt]
        return res[:conv] - 1.0
    end
    alpha_opt = secant(root_fun, alpha, ustart)

    sol = optimizer(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol
end

"""
    bryan
"""
function bryan(mec::MaxEntContext)
    println("Using bryan algorithm to solve the maximum entropy problem")

    use_bayes = true
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    maxprob = 0.0
    while true
        sol = optimizer(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / ratio
        @. ustart = sol[:u_opt]
        prob = sol[:prob]
        if prob > maxprob
            maxprob = prob
        elseif prob < 0.01 * maxprob
            break
        end
    end

    alpharr = map(x->x[:alpha], optarr)
    probarr = map(x->x[:prob], optarr)
    specarr = map(x->x[:A_opt], optarr)
    probarr = probarr ./ (-trapz(alpharr, probarr))

    nprob = length(probarr)
    nmesh = length(specarr[1])
    A_opt = zeros(F64, nmesh)
    Spectra = zeros(F64, nprob, nmesh)
    for i = 1:nprob
        Spectra[i,:] = specarr[i] * probarr[i]
    end
    for i = 1:nmesh
        A_opt[i] = -trapz(alpharr, Spectra[:,i])
    end

    sol = Dict{Symbol,Any}()
    sol[:A_opt] = A_opt

    return optarr, sol
end

"""
    chi2kink
"""
function chi2kink(mec::MaxEntContext)
    function fitfun(x, p)
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    println("Apply chi2kink algorithm to determine optimized α")

    use_bayes = false
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    nalph = get_m("nalph")
    alpha_end = alpha / (ratio^nalph)
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    chi_arr = []
    alpharr = []
    while true
        sol = optimizer(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        push!(chi_arr, sol[:chi2])
        push!(alpharr, alpha)
        @. ustart = sol[:u_opt]
        alpha = alpha / ratio
        if alpha < alpha_end
            break
        end
    end

    good_data = isfinite.(chi_arr)
    fit = curve_fit(fitfun, log10.(alpharr[good_data]), log10.(chi_arr[good_data]), [0.0, 5.0, 2.0, 0.0])
    a, b, c, d = fit.param

    fit_position = 2.5
    a_opt = c - fit_position / d
    closest_idx = argmin( abs.( log10.(alpharr) .- a_opt ) )
    ustart = optarr[closest_idx][:u_opt]
    alpha_opt = 10.0 ^ a_opt

    sol = optimizer(mec, alpha_opt, ustart, use_bayes)

    open("chi2kink.fit", "w") do fout
        alpharr = log10.(alpharr)
        chi_arr = log10.(chi_arr)
        fit_arr = fitfun(alpharr, fit.param)
        for i = 1:length(alpharr)
            println(fout, alpharr[i], " ", chi_arr[i], " ", fit_arr[i])
        end
    end

    return optarr, sol
end

"""
    optimizer
"""
function optimizer(mec::MaxEntContext, alpha::F64, ustart::Vector{F64}, use_bayes::Bool)
    blur = get_m("blur")
    offdiag = get_c("offdiag")

    if offdiag
        solution, call = newton(f_and_J_offdiag, ustart, mec, alpha)
        u_opt = copy(solution)
        A_opt = svd_to_real_offdiag(mec, solution)
        entropy = calc_entropy_offdiag(mec, A_opt, u_opt)
    else
        solution, call = newton(f_and_J, ustart, mec, alpha)
        u_opt = copy(solution)
        A_opt = svd_to_real(mec, solution)
        entropy = calc_entropy(mec, A_opt, u_opt)
    end

    χ² = calc_chi2(mec, A_opt)
    norm = trapz(mec.mesh, A_opt)

    result_dict = Dict{Symbol,Any}()
    result_dict[:u_opt] = u_opt
    result_dict[:alpha] = alpha
    result_dict[:entropy] = entropy
    result_dict[:chi2] = χ²
    result_dict[:norm] = norm
    result_dict[:Q] = alpha * entropy - 0.5 * χ²
    if blur > 0.0
        make_blur(mec.mesh, A_opt, blur)
        result_dict[:A_opt] = A_opt
    else
        result_dict[:A_opt] = A_opt
    end

    if use_bayes
        if offdiag
            ng, tr, conv, prob = calc_bayes_offdiag(mec, A_opt, entr, χ², alpha)
        else
            ng, tr, conv, prob = calc_bayes(mec, A_opt, entr, χ², alpha)
        end
        result_dict[:n_good] = ng
        result_dict[:trace] = tr
        result_dict[:conv] = conv
        result_dict[:prob] = prob
    end

    @printf("log10(α) = %8.4f ", log10(alpha))
    @printf("χ² = %8.4e ", χ²)
    @printf("S = %8.4e ", entropy)
    @printf("call = %4i ", call)
    @printf("norm = %8.4f\n", norm)

    return result_dict
end

"""
    precompute
"""
function precompute(Gᵥ::Vector{F64}, σ²::Vector{F64}, mesh::AbstractMesh, model::Vector{F64}, kernel::Matrix{F64})
    offdiag = get_c("offdiag")

    U, V, S = make_singular_space(kernel)

    nmesh = mesh.nmesh
    weight = mesh.weight
    n_svd = length(S)

    W₂ = zeros(F64, n_svd, nmesh)
    W₃ = zeros(F64, n_svd, n_svd, nmesh)
    Bₘ = zeros(F64, n_svd)
    hess = zeros(F64, nmesh, nmesh)

    if offdiag
        @einsum W₂[m,l] = σ²[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * weight[l]
    else
        @einsum W₂[m,l] = σ²[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * weight[l] * model[l]
    end

    for i = 1:nmesh
        W₃[:,:,i] = W₂[:,i] * adjoint(V[i,:])
    end

    @einsum Bₘ[m] = S[m] * U[k,m] * σ²[k] * Gᵥ[k]

    @einsum hess[i,j] = weight[i] * weight[j] * kernel[k,i] * kernel[k,j] * σ²[k]

    return V, W₂, W₃, Bₘ, hess
end

"""
    f_and_J
"""
function f_and_J(u::Vector{F64}, mec::MaxEntContext, alpha::F64)
    v = mec.Vₛ * u
    w = exp.(v)

    n_svd = length(mec.Bₘ)
    term_1 = zeros(F64, n_svd)
    term_2 = zeros(F64, n_svd, n_svd)

    for i = 1:n_svd
        term_1[i] = dot(mec.W₂[i,:], w)
    end

    for j = 1:n_svd
        for i = 1:n_svd
            term_2[i,j] = dot(mec.W₃[i,j,:], w)
        end
    end

    f = alpha * u + term_1 - mec.Bₘ
    J = alpha * diagm(ones(n_svd)) + term_2

    return f, J
end

"""
    f_and_J_offdiag
"""
function f_and_J_offdiag(u::Vector{F64}, mec::MaxEntContext, alpha::F64)
    v = mec.Vₛ * u
    w = exp.(v)

    a_plus = mec.model .* w
    a_minus = mec.model ./ w
    a1 = a_plus - a_minus
    a2 = a_plus + a_minus

    n_svd = length(mec.Bₘ)
    term_1 = zeros(F64, n_svd)
    term_2 = zeros(F64, n_svd, n_svd)

    for i = 1:n_svd
        term_1[i] = dot(mec.W₂[i,:], a1)
    end

    for j = 1:n_svd
        for i = 1:n_svd
            term_2[i,j] = dot(mec.W₃[i,j,:], a2)
        end
    end

    f = alpha * u + term_1 - mec.Bₘ
    J = alpha * diagm(ones(n_svd)) + term_2

    return f, J
end

"""
    svd_to_real
"""
function svd_to_real(mec::MaxEntContext, u::Vector{F64})
    return mec.model .* exp.(mec.Vₛ * u)
end

"""
    svd_to_real_offdiag
"""
function svd_to_real_offdiag(mec::MaxEntContext, u::Vector{F64})
    w = exp.(mec.Vₛ * u)
    return (mec.model .* w) - (mec.model ./ w)
end

"""
    calc_entropy
"""
function calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})
    f = A - mec.model - A .* (mec.Vₛ * u)
    return trapz(mec.mesh, f)
end

"""
    calc_entropy_offdiag
"""
function calc_entropy_offdiag(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})
    root = sqrt.(A .^ 2.0 + 4.0 .* mec.model .* mec.model)
    f = root - mec.model - mec.model - A .* log.((root + A) ./ (2.0 * mec.model))
    return trapz(mec.mesh, f)
end

"""
    calc_bayes
"""
function calc_bayes(mec::MaxEntContext, A::Vector{F64}, ent::F64, chisq::F64, alpha::F64)
    mesh = mec.mesh

    T = sqrt.(A ./ mesh.weight)
    Λ = (T * T') .* mec.hess

    lam = eigvals(Hermitian(Λ))
    ng = -2.0 * alpha * ent
    tr = sum(lam ./ (alpha .+ lam))
    conv = tr / ng

    eig_sum = sum(log.(alpha ./ (alpha .+ lam)))
    log_prob = alpha * ent - 0.5 * chisq + log(alpha) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

"""
    calc_bayes_offdiag
"""
function calc_bayes_offdiag(mec::MaxEntContext, A::Vector{F64}, ent::F64, chisq::F64, alpha::F64)
    mesh = mec.mesh

    T = (( A .^ 2.0 + 4.0 * mec.model .* mec.model ) / (mesh.weight .^ 2.0)) .^ 0.25
    Λ = (T * T') .* mec.hess

    lam = eigvals(Hermitian(Λ))
    ng = -2.0 * alpha * ent
    tr = sum(lam ./ (alpha .+ lam))
    conv = tr / ng

    eig_sum = sum(log.(alpha ./ (alpha .+ lam)))
    log_prob = alpha * ent - 0.5 * chisq + log(alpha) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

"""
    calc_chi2
"""
function calc_chi2(mec::MaxEntContext, A::Vector{F64})
    ndim, _ = size(mec.kernel)

    T = zeros(F64, ndim)
    for i = 1:ndim
        T[i] = mec.Gᵥ[i] - trapz(mec.mesh, mec.kernel[i,:] .* A)
    end
    value = sum(mec.σ² .* T .^ 2.0)

    return value
end
