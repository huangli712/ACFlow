#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/10/10
#

#=
### *Customized Structs* : *NevanAC Solver*
=#

"""
    NevanACContext

Mutable struct. It is used within the NevanAC solver only.

### Members

* Gᵥ   -> Input data for correlator.
* grid -> Grid for input data.
* mesh -> Mesh for output spectrum.
* Φ    -> `Φ` vector in Schur algorithm.
* 𝒜    -> Coefficients matrix `abcd` in Schur algorithm.
* ℋ    -> Hardy matrix for Hardy basis optimization.
* 𝑎𝑏   -> Coefficients matrix for expanding `θ` with Hardy basis.
* hmin -> Minimal value of the order of Hardy basis functions.
* hopt -> Optimal value of the order of Hardy basis functions.
"""
mutable struct NevanACContext
    Gᵥ   :: Vector{APC}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    Φ    :: Vector{APC}
    𝒜    :: Array{APC,3}
    ℋ    :: Array{APC,2}
    𝑎𝑏   :: Vector{C64}
    hmin :: I64
    hopt :: I64
end

#=
### *Global Drivers*
=#

"""
    solve(S::NevanACSolver, rd::RawData)

Solve the analytic continuation problem by the Nevanlinna analytical
continuation method.
"""
function solve(S::NevanACSolver, rd::RawData)
    println("[ NevanAC ]")
    nac = init(S, rd)
    run(nac)
    Aout, Gout = last(nac)
    return nac.mesh.mesh, Aout, Gout
end

"""
    init(S::NevanACSolver, rd::RawData)

Initialize the NevanAC solver and return a NevanACContext struct.
"""
function init(S::NevanACSolver, rd::RawData)
    # Setup numerical precision. Note that the NAC method is extremely
    # sensitive to the float point precision.
    setprecision(128)

    # Convert the input data to APC, i.e., Complex{BigFloat}.
    ωₙ = APC.(rd._grid * im)
    Gₙ = APC.(rd.value)

    # Evaluate the optimal value for the size of input data.
    # Here we just apply the Pick criterion.
    ngrid = calc_noptim(ωₙ, Gₙ)

    # Prepera input data
    Gᵥ = calc_mobius(-Gₙ[1:ngrid])
    reverse!(Gᵥ)
    println("Postprocess input data: ", length(Gᵥ), " points")

    # Prepare grid for input data
    grid = make_grid(rd, T = APF)
    resize!(grid, ngrid)
    reverse!(grid)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh(T = APF)
    println("Build mesh for spectrum: ", length(mesh), " points")

    # Precompute key quantities to accelerate the computation
    Φ, 𝒜, ℋ, 𝑎𝑏 = precompute(grid, mesh, Gᵥ)
    println("Precompute key matrices")

    # Actually, the NevanACContext struct already contains enough
    # information to build the Nevanlinna interpolant and get the
    # spectrum, but Hardy basis optimization is needed to smooth
    # the results further.
    return NevanACContext(Gᵥ, grid, mesh, Φ, 𝒜, ℋ, 𝑎𝑏, 1, 1)
end

"""
    run(nac::NevanACContext)

Perform Hardy basis optimization to smooth the spectrum. the members `ℋ`,
`𝑎𝑏`, `hmin`, and `hopt` of the NevanACContext object (`nac`) should be
updated in this function.
"""
function run(nac::NevanACContext)
    hardy = get_n("hardy")
    #
    if hardy
        # Determine the minimal Hardy order (`hmin`), update `ℋ` and `𝑎𝑏`.
        calc_hmin!(nac)

        # Determine the optimal Hardy order (`hopt`), update `ℋ` and `𝑎𝑏`.
        calc_hopt!(nac)
    end
end

"""
    last(nac::NevanACContext)

Postprocess the results generated during the Nevanlinna analytical
continuation simulations.
"""
function last(nac::NevanACContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Calculate full response function on real axis and write them
    # Note that _G is actually 𝑁G, so there is a `-` symbol for the
    # return value.
    _G = C64.(calc_green(nac.𝒜, nac.ℋ, nac.𝑎𝑏))
    fwrite && write_complete(nac.mesh, -_G)

    # Calculate and write the spectral function
    Aout = F64.(imag.(_G) ./ π)
    fwrite && write_spectrum(nac.mesh, Aout)

    # Regenerate the input data and write them
    kernel = make_kernel(nac.mesh, nac.grid)
    G = reprod(nac.mesh, kernel, Aout)
    fwrite && write_backward(nac.grid, G)

    return Aout, -_G
end

#=
### *Service Functions*
=#

#=
*Remarks* :

**Mobius transformation**

```math
\begin{equation}
z \mapsto \frac{z - i}{z + i}.
\end{equation}
```

It maps `z` from ``\overline{\mathcal{C}^{+}}`` to ``\overline{\mathcal{D}}``.
See `calc_mobius()`.

---

**Inverse Mobius transformation**

```math
\begin{equation}
z \mapsto i \frac{1 + z}{1 - z}.
\end{equation}
```

It maps `z` from ``\overline{\mathcal{D}}`` to ``\overline{\mathcal{C}^{+}}``.
See `calc_inv_mobius()`.

---

**Pick matrix**

```math
\begin{equation}
\mathcal{P} = 
\left[
    \frac{1-\lambda_i \lambda^*_j}{1-h(Y_i)h(Y_j)^*}
\right]_{i,j}.
\end{equation}
```

Here, ``Y_i`` is the *i*th Matsubara frequency, ``C_i`` is the value of
``\mathcal{NG}`` at ``Y_i``, and ``\lambda_i`` is the value of ``\theta``
at ``Y_i``:

```math
\begin{equation}
\lambda_i = \theta(Y_i) = \frac{C_i - i}{C_i + i},
\end{equation}
```

```math
\begin{equation}
h(Y_i) = \frac{Y_i - i}{Y_i + i}.
\end{equation}
```

See `calc_pick()`.

---

**Φ vector**

At first, we have

```math
\begin{equation}
\phi_{\alpha} = \theta_{\alpha}(Y_{\alpha}).
\end{equation}
```

So

```math
\begin{equation}
\phi_1 = \theta_1 (Y_1),
\end{equation}
```

```math
\begin{equation}
\phi_{\beta} =
\frac{-d_{\beta}\theta(Y_{\beta}) + b_{\beta}}{c_{\beta}\theta(Y_{\beta}) - \alpha_{\beta}}
\end{equation}
```

where

```math
\begin{equation}
\begin{pmatrix}
a_{\beta} & b_{\beta} \\
c_{\beta} & d_{\beta}
\end{pmatrix}
= \prod^{\beta-1}_{\alpha=1}
\begin{pmatrix}
\frac{Y_{\beta}-Y_{\alpha}}{Y_{\beta}-Y^*_{\alpha}} & \phi_{\alpha} \\
\phi^*_{\alpha}\frac{Y_{\beta} - Y_{\alpha}}{Y_{\beta} - Y^*_{\alpha}} & 1
\end{pmatrix},
\end{equation}
```

See `calc_phis()`.

---

**Contractive function θ(z)**

```math
\begin{equation}
\theta(z)[z; θ_{M+1}(z)] =
\frac{a(z)\theta_{M+1}(z) + b(z)}{c(z)\theta_{M+1}(z) + d(z)}
\end{equation}
```

where

```math
\begin{equation}
\begin{pmatrix}
a(z) & b(z) \\
c(z) & d(z) 
\end{pmatrix}
= \prod^{M}_{j=1}
\begin{pmatrix}
\frac{z - Y_j}{z - Y^*_j} & \phi_j \\
\phi^*_j \frac{z - Y_j}{z - Y^*_j} & 1
\end{pmatrix}
\end{equation}
```

See `calc_abcd()` and `calc_theta()`.

---

**Hardy basis**

```math
\begin{equation}
f^k(z) = \frac{1}{\sqrt{\pi}(z + i)}
    \left( \frac{z - i}{z + i} \right)^k.
\end{equation}
```

See `calc_hbasis()` and `calc_hmatrix()`.

---

**Expanding ``\theta_{M+1}`` using Hardy basis**

```math
\theta_{M+1} = \sum^{H}_{k=0} \left[a_k f^k(z) + b_k f^k(z)^*\right]
```

See `calc_theta()`.

---

**Smooth norm**

```math
\begin{equation}
F[A_{\theta_{M+1}}(\omega)] =
    \left|
        1 - \int A_{\theta_{M+1}}(\omega) d\omega
    \right|^2 +
    \alpha \int \left[A^{''}_{\theta_{M+1}}(\omega)\right]^2 d\omega,
\end{equation}
```
where the first term enforces proper normalization while the second
term promotes smoothness by minimizing second derivatives. Here α
is a regulation parameter.

See `smooth_norm()`.
=#

"""
    precompute(grid::AbstractGrid,
               mesh::AbstractMesh,
               Gᵥ::Vector{APC})

Precompute some key quantities, such as `Φ`, `𝒜`, `ℋ`, and `𝑎𝑏`. Note
that `Φ` and `𝒜` won't be changed any more. But `𝒜` and `𝑎𝑏` should be
updated by the Hardy basis optimization to get a smooth spectrum. Here
`Gᵥ` is input data, `grid` is the grid for input data, and `mesh` is
the mesh for output spectrum.
"""
function precompute(grid::AbstractGrid,
                    mesh::AbstractMesh,
                    Gᵥ::Vector{APC})
    # Evaluate ϕ and `abcd` matrices
    Φ = calc_phis(grid, Gᵥ)
    𝒜 = calc_abcd(grid, mesh, Φ)

    # Allocate memory for evaluating θ
    # The initial Hardy order is just 1.
    ℋ = calc_hmatrix(mesh, 1)
    𝑎𝑏 = zeros(C64, 2)

    return Φ, 𝒜, ℋ, 𝑎𝑏
end

"""
    calc_mobius(z::Vector{APC})

A direct Mobius transformation.
"""
function calc_mobius(z::Vector{APC})
    return @. (z - im) / (z + im)
end

"""
    calc_inv_mobius(z::Vector{APC})

An inverse Mobius transformation.
"""
function calc_inv_mobius(z::Vector{APC})
    return @. im * (one(APC) + z) / (one(APC) - z)
end

"""
    calc_pick(k::I64, ℎ::Vector{APC}, λ::Vector{APC})

Try to calculate the Pick matrix, anc check whether it is a positive
semidefinite matrix. See Eq. (5) in Fei's NAC paper.
"""
function calc_pick(k::I64, ℎ::Vector{APC}, λ::Vector{APC})
    pick = zeros(APC, k, k)

    # Calculate the Pick matrix
    for j = 1:k
        for i = 1:k
            num = one(APC) - λ[i] * conj(λ[j])
            den = one(APC) - ℎ[i] * conj(ℎ[j])
            pick[i,j] = num / den
        end
        pick[j,j] += APC(1e-250)
    end

    # Cholesky decomposition
    return issuccess(cholesky(pick, check = false))
end

"""
    calc_phis(grid::AbstractGrid, Gᵥ::Vector{APC})

Try to calculate the Φ vector, which is used to calculate the 𝒜 matrix.
Note that Φ should not be changed anymore once it has been established.
"""
function calc_phis(grid::AbstractGrid, Gᵥ::Vector{APC})
    ngrid = length(grid)

    # Allocate memory
    Φ = zeros(APC, ngrid) 
    𝒜 = zeros(APC, 2, 2, ngrid)
    ∏ = zeros(APC, 2, 2)
    𝑔 = grid.ω * im

    # Initialize the `abcd` matrix
    for i = 1:ngrid
        𝒜[:,:,i] .= Matrix{APC}(I, 2, 2)
    end

    # Evaluate Φ using recursive algorithm
    Φ[1] = Gᵥ[1]
    for j = 1:ngrid-1
        for k = j+1:ngrid
            ∏[1,1] = ( 𝑔[k] - 𝑔[j] ) / ( 𝑔[k] - conj(𝑔[j]) )
            ∏[1,2] = Φ[j]
            ∏[2,1] = conj(Φ[j]) * ∏[1,1]
            ∏[2,2] = one(APC)
            view(𝒜,:,:,k) .= view(𝒜,:,:,k) * ∏
        end
        num = 𝒜[1,2,j+1] - 𝒜[2,2,j+1] * Gᵥ[j+1]
        den = 𝒜[2,1,j+1] * Gᵥ[j+1] - 𝒜[1,1,j+1]
        Φ[j+1] = num / den
    end

    return Φ
end

"""
    calc_abcd(grid::AbstractGrid, mesh::AbstractMesh, Φ::Vector{APC})

Try to calculate the coefficients matrix abcd (here it is called 𝒜),
which is then used to calculate θ. See Eq. (8) in Fei's NAC paper.
"""
function calc_abcd(grid::AbstractGrid, mesh::AbstractMesh, Φ::Vector{APC})
    eta::APF = get_n("eta")

    ngrid = length(grid)
    nmesh = length(mesh)

    𝑔 = grid.ω * im
    𝑚 = mesh.mesh .+ im * eta

    𝒜 = zeros(APC, 2, 2, nmesh)
    ∏ = zeros(APC, 2, 2)

    for i = 1:nmesh
        result = Matrix{APC}(I, 2, 2)
        𝑧 = 𝑚[i]
        for j = 1:ngrid
            ∏[1,1] = ( 𝑧 - 𝑔[j] ) / ( 𝑧 - conj(𝑔[j]) )
            ∏[1,2] = Φ[j]
            ∏[2,1] = conj(Φ[j]) * ∏[1,1]
            ∏[2,2] = one(APC)
            result *= ∏
        end

        𝒜[:,:,i] .= result
    end

    return 𝒜
end

"""
    calc_hbasis(z::APC, k::I64)

Try to calculate the Hardy basis ``f^k(z)``.
"""
function calc_hbasis(z::APC, k::I64)
    w = ( z - im ) / ( z + im )
    return 1.0 / ( sqrt(π) * (z + im) ) * w^k
end

"""
    calc_hmatrix(mesh::AbstractMesh, H::I64)

Try to calculate ``[f^k(z), f^k(z)^*]`` for 0 ≤ 𝑘 ≤ 𝐻-1, which is
called the hardy matrix (ℋ) and is used to evaluate ``\theta_{M+1}``.
"""
function calc_hmatrix(mesh::AbstractMesh, H::I64)
    # Build real axis
    eta::APF = get_n("eta")
    𝑚 = mesh.mesh .+ eta * im

    # Allocate memory for the Hardy matrix
    nmesh = length(mesh)
    ℋ = zeros(APC, nmesh, 2*H)

    # Build the Hardy matrix
    for k = 1:H
        ℋ[:,2*k-1] .= calc_hbasis.(𝑚,k-1)
        ℋ[:,2*k]   .= conj(ℋ[:,2*k-1])
    end

    return ℋ
end

"""
    calc_theta(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})

Try to calculate the contractive function θ(z). 𝒜 is the coefficients
matrix abcd, ℋ is the Hardy matrix, and 𝑎𝑏 are complex coefficients
for expanding θₘ₊₁. See Eq. (7) in Fei's NAC paper.
"""
function calc_theta(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
    # Well, we should calculate θₘ₊₁ at first.
    θₘ₊₁ = ℋ * 𝑎𝑏

    # Then we evaluate θ according Eq. (7)
    num = 𝒜[1,1,:] .* θₘ₊₁ .+ 𝒜[1,2,:]
    den = 𝒜[2,1,:] .* θₘ₊₁ .+ 𝒜[2,2,:]
    θ = num ./ den

    return θ
end

"""
    calc_green(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})

Firstly we try to calculate θ. Then θ is back transformed to a Nevanlinna
interpolant via the inverse Mobius transform. Here, `𝒜` (`abcd` matrix),
`ℋ` (Hardy matrix), and `𝑎𝑏` are used to evaluate θ.
"""
function calc_green(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
    θ = calc_theta(𝒜, ℋ, 𝑎𝑏)
    gout = calc_inv_mobius(θ)

    return gout
end

"""
    calc_noptim(ωₙ::Vector{APC}, Gₙ::Vector{APC})

Evaluate the optimal value for the size of input data (how may frequency
points are actually used in the analytic continuation simulations) via
the Pick criterion.
"""
function calc_noptim(ωₙ::Vector{APC}, Gₙ::Vector{APC})
    # Get size of input data
    ngrid = length(ωₙ)

    # Check whether the Pick criterion is applied 
    pick = get_n("pick")
    if !pick
        return ngrid
    end

    # Apply invertible Mobius transformation. We actually work at
    # the \bar{𝒟} space.
    𝓏 = calc_mobius(ωₙ)
    𝒢 = calc_mobius(-Gₙ)

    # Find the optimal value of k until the Pick criterion is violated
    k = 0
    success = true
    while success && k ≤ ngrid
        k += 1
        success = calc_pick(k, 𝓏, 𝒢)
    end

    # Return the optimal value for the size of input data
    if !success
        println("The size of input data is optimized to $(k-1)")
        return k - 1
    else
        println("The size of input data is optimized to $(ngrid)")
        return ngrid
    end
end

"""
    calc_hmin!(nac::NevanACContext)

Try to perform Hardy basis optimization. Such that the Hardy matrix ℋ
and the corresponding coefficients 𝑎𝑏 are updated. They are used to
calculate θ, which is then back transformed to generate smooth G (i.e.,
the spectrum) at real axis.

This function will determine the minimal value of H (hmin). Of course,
ℋ and 𝑎𝑏 in NevanACContext object are also changed.
"""
function calc_hmin!(nac::NevanACContext)
    hmax = get_n("hmax")

    h = 1
    while h ≤ hmax
        println("H = $h")

        # Prepare initial ℋ and 𝑎𝑏
        ℋ = calc_hmatrix(nac.mesh, h)
        𝑎𝑏 = zeros(C64, 2*h)

        # Hardy basis optimization
        causality, optim = hardy_optimize!(nac, ℋ, 𝑎𝑏, h)

        # Check whether the causality is preserved and the
        # optimization is successful.
        if causality && optim
            nac.hmin = h
            break
        else
            h = h + 1
        end
    end
end

"""
    calc_hopt!(nac::NevanACContext)

Try to perform Hardy basis optimization. Such that the Hardy matrix ℋ
and the corresponding coefficients 𝑎𝑏 are updated. They are used to
calculate θ, which is then back transformed to generate smooth G (i.e.,
the spectrum) at real axis.

This function will determine the optimal value of H (hopt). Of course,
ℋ and 𝑎𝑏 in NevanACContext object are also changed.
"""
function calc_hopt!(nac::NevanACContext)
    hmax = get_n("hmax")
    hmax = 3 # DEBUG

    for h = nac.hmin + 1:hmax
        println("H = $h")

        # Prepare initial ℋ and 𝑎𝑏
        ℋ = calc_hmatrix(nac.mesh, h)
        𝑎𝑏  = copy(nac.𝑎𝑏)
        push!(𝑎𝑏, zero(C64))
        push!(𝑎𝑏, zero(C64))
        @assert size(ℋ)[2] == length(𝑎𝑏)

        # Hardy basis optimization
        causality, optim = hardy_optimize!(nac, ℋ, 𝑎𝑏, h)

        # Check whether the causality is preserved and the
        # optimization is successful.
        if !(causality && optim)
            break
        end
    end
end

"""
    hardy_optimize!(nac::NevanACContext,
                    ℋ::Array{APC,2},
                    𝑎𝑏::Vector{C64},
                    H::I64)

For given Hardy matrix ℋ, try to update the expanding coefficients 𝑎𝑏
by minimizing the smooth norm.
"""
function hardy_optimize!(nac::NevanACContext,
                         ℋ::Array{APC,2},
                         𝑎𝑏::Vector{C64},
                         H::I64)
    function 𝑓(x::Vector{C64})
        return smooth_norm(nac, ℋ, x)
    end

    function 𝐽!(J::Vector{C64}, x::Vector{C64})
        J .= gradient(𝑓, x)
    end

    res = optimize(𝑓, 𝐽!, 𝑎𝑏, BFGS(), 
                   Options(iterations = 500, show_trace = true))
    
    @show res.minimizer

    if  !(converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(ℋ, res.minimizer)

    if causality && (converged(res))
        nac.hopt = H
        nac.𝑎𝑏 = res.minimizer
        nac.ℋ = ℋ
    end
    
    return causality, (converged(res))
end

"""
    smooth_norm(nac::NevanACContext, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})

Establish the smooth norm, which is used to improve the smoothness of
the output spectrum.
"""
function smooth_norm(nac::NevanACContext, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
    # Get regulation parameter
    α = get_n("alpha")

    # Generate output spectrum
    _G = calc_green(nac.𝒜, ℋ, 𝑎𝑏)
    A = F64.(imag.(_G) ./ π)

    # Normalization term
    𝑓₁ = trapz(nac.mesh, A)

    # Smoothness term
    sd = deriv2(nac.mesh.mesh, A)
    x_sd = nac.mesh.mesh[2:end-1]
    𝑓₂ = trapz(x_sd, abs.(sd) .^ 2)

    # Assemble the final smooth norm
    𝐹 = abs(1.0 - 𝑓₁)^2 + α * 𝑓₂

    return F64(𝐹)
end

"""
"""
function check_pick(wn::Vector{APC}, gw::Vector{APC}, Nopt::I64)
    freq = calc_mobius(wn[1:Nopt])
    val = calc_mobius(-gw[1:Nopt])

    success = calc_pick(Nopt, val, freq)
    if success
        println("Pick matrix is positive semi-definite.")
    else
        println("Pick matrix is non positive semi-definite matrix in Schur method.")
    end
    
    freq = reverse(wn[1:Nopt])
    val  = reverse(val)
end

"""
"""
function check_causality(ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
    param = ℋ * 𝑎𝑏
    @show typeof(param), size(ℋ), size(𝑎𝑏), size(param)

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        println("max_theta = ",max_theta)
        println("Hardy optimization was success.")
        causality = true
    else
        println("max_theta = ",max_theta)
        println("Hardy optimization was failure.")
        causality = false
    end

    return causality
end

struct Manifold end
abstract type AbstractOptimizer end

struct OptimizationState{Tf<:Real}
    iteration::Int
    value::Tf
    g_norm::Tf
    metadata::Dict
end

struct BFGS{IL, L, H, T, TM} <: AbstractOptimizer
    alphaguess!::IL
    linesearch!::L
    initial_invH::H
    initial_stepnorm::T
    manifold::TM
end

mutable struct BFGSState{Tx, Tm, T,G}
    x::Tx
    x_previous::Tx
    g_previous::G
    f_x_previous::T
    dx::Tx
    dg::Tx
    u::Tx
    invH::Tm
    s::Tx
    x_ls::Tx
    alpha::T
end

abstract type AbstractObjective end
mutable struct ManifoldObjective{T <: AbstractObjective} <: AbstractObjective
    manifold::Manifold
    inner_obj::T
end
# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable1{TF, TDF, TX} <: AbstractObjective
    f # objective
    df # (partial) derivative of objective
    fdf # objective and (partial) derivative of objective
    F::TF # cache for f output
    DF::TDF # cache for df output
    x_f::TX # x used to evaluate f (stored in F)
    x_df::TX # x used to evaluate df (stored in DF)
    f_calls::Vector{Int}
    df_calls::Vector{Int}
end

mutable struct MultivariateOptimizationResults{O, Tx, Tc, Tf, M, Tls, Tsb}
    method::O
    initial_x::Tx
    minimizer::Tx
    minimum::Tf
    iterations::Int
    iteration_converged::Bool
    x_converged::Bool
    x_abstol::Tf
    x_reltol::Tf
    x_abschange::Tc
    x_relchange::Tc
    f_converged::Bool
    f_abstol::Tf
    f_reltol::Tf
    f_abschange::Tc
    f_relchange::Tc
    g_converged::Bool
    g_abstol::Tf
    g_residual::Tc
    f_increased::Bool
    trace::M
    f_calls::Int
    g_calls::Int
    h_calls::Int
    ls_success::Tls
    time_limit::Float64
    time_run::Float64
    stopped_by::Tsb
end

struct Options{T, TCallback}
    x_abstol::T
    x_reltol::T
    f_abstol::T
    f_reltol::T
    g_abstol::T
    g_reltol::T
    outer_x_abstol::T
    outer_x_reltol::T
    outer_f_abstol::T
    outer_f_reltol::T
    outer_g_abstol::T
    outer_g_reltol::T
    f_calls_limit::Int
    g_calls_limit::Int
    h_calls_limit::Int
    allow_f_increases::Bool
    allow_outer_f_increases::Bool
    successive_f_tol::Int
    iterations::Int
    outer_iterations::Int
    store_trace::Bool
    #trace_simplex::Bool
    show_trace::Bool
    extended_trace::Bool
    show_every::Int
    callback::TCallback
    time_limit::Float64
end

function Options(;
        x_tol = nothing,
        f_tol = nothing,
        g_tol = nothing,
        x_abstol::Real = 0.0,
        x_reltol::Real = 0.0,
        f_abstol::Real = 0.0,
        f_reltol::Real = 0.0,
        g_abstol::Real = 1e-8,
        g_reltol::Real = 1e-8,
        outer_x_tol = nothing,
        outer_f_tol = nothing,
        outer_g_tol = nothing,
        outer_x_abstol::Real = 0.0,
        outer_x_reltol::Real = 0.0,
        outer_f_abstol::Real = 0.0,
        outer_f_reltol::Real = 0.0,
        outer_g_abstol::Real = 1e-8,
        outer_g_reltol::Real = 1e-8,
        f_calls_limit::Int = 0,
        g_calls_limit::Int = 0,
        h_calls_limit::Int = 0,
        allow_f_increases::Bool = true,
        allow_outer_f_increases::Bool = true,
        successive_f_tol::Int = 1,
        iterations::Int = 1_000,
        outer_iterations::Int = 1000,
        store_trace::Bool = false,
        #trace_simplex::Bool = false,
        show_trace::Bool = false,
        extended_trace::Bool = false,
        show_every::Int = 1,
        callback = nothing,
        time_limit = NaN)
    show_every = show_every > 0 ? show_every : 1
    #if extended_trace && callback === nothing
    #    show_trace = true
    #end
    if !(x_tol === nothing)
        x_abstol = x_tol
    end
    if !(g_tol === nothing)
        g_abstol = g_tol
    end
    if !(f_tol === nothing)
        f_reltol = f_tol
    end
    if !(outer_x_tol === nothing)
        outer_x_abstol = outer_x_tol
    end
    if !(outer_g_tol === nothing)
        outer_g_abstol = outer_g_tol
    end
    if !(outer_f_tol === nothing)
        outer_f_reltol = outer_f_tol
    end
    Options(promote(x_abstol, x_reltol, f_abstol, f_reltol, g_abstol, g_reltol, outer_x_abstol, outer_x_reltol, outer_f_abstol, outer_f_reltol, outer_g_abstol, outer_g_reltol)..., f_calls_limit, g_calls_limit, h_calls_limit,
        allow_f_increases, allow_outer_f_increases, successive_f_tol, Int(iterations), Int(outer_iterations), store_trace, show_trace, extended_trace,
        Int(show_every), callback, Float64(time_limit))
end

const OptimizationTrace{Tf} = Vector{OptimizationState{Tf}}

include("hagerzhang.jl")

function BFGS(; alphaguess = InitialStatic(),
    linesearch = HagerZhang(),
    initial_invH = nothing,
    initial_stepnorm = nothing,
    manifold::Manifold=Manifold())
    BFGS(alphaguess, linesearch, initial_invH, initial_stepnorm, manifold)
end

function OnceDifferentiable1(f, df,
                   x::AbstractArray,
                   F::Real = real(zero(eltype(x))),
                   DF::AbstractArray = alloc_DF(x, F))
    function fdf(gx, x)
        df(gx, x)
        return f(x)
    end
    x_f, x_df = x_of_nans(x), x_of_nans(x)
    OnceDifferentiable1(f, df, fdf, copy(F), copy(DF), x_f, x_df, [0,], [0,])
end

function initial_state(method::BFGS, d, initial_x::AbstractArray{T}) where T
    initial_x = copy(initial_x)
    retract!(method.manifold, initial_x)

    value_gradient!!(d, initial_x)

    project_tangent!(method.manifold, gradient(d), initial_x)

    if method.initial_invH === nothing
        if method.initial_stepnorm === nothing
            # Identity matrix of size n x n
            invH0 = _init_identity_matrix(initial_x)
        else
            initial_scale = T(method.initial_stepnorm) * inv(norm(gradient(d), Inf))
            invH0 = _init_identity_matrix(initial_x, initial_scale)
        end
    else
        invH0 = method.initial_invH(initial_x)
    end
    # Maintain a cache for line search results
    # Trace the history of states visited
    BFGSState(initial_x, # Maintain current state in state.x
              copy(initial_x), # Maintain previous state in state.x_previous
              copy(gradient(d)), # Store previous gradient in state.g_previous
              real(T)(NaN), # Store previous f in state.f_x_previous
              similar(initial_x), # Store changes in position in state.dx
              similar(initial_x), # Store changes in gradient in state.dg
              similar(initial_x), # Buffer stored in state.u
              invH0, # Store current invH in state.invH
              similar(initial_x), # Store current search direction in state.s
              similar(initial_x), # Buffer of x for line search in state.x_ls
              real(one(T))
            )
end

function update_state!(d, state::BFGSState, method::BFGS)
    n = length(state.x)
    T = eltype(state.s)
    # Set the search direction
    # Search direction is the negative gradient divided by the approximate Hessian
    mul!(vec(state.s), state.invH, vec(gradient(d)))
    rmul!(state.s, T(-1))
    project_tangent!(method.manifold, state.s, state.x)

    # Maintain a record of the previous gradient
    copyto!(state.g_previous, gradient(d))

    # Determine the distance of movement along the search line
    # This call resets invH to initial_invH is the former in not positive
    # semi-definite
    lssuccess = perform_linesearch!(state, method, ManifoldObjective(method.manifold, d))

    # Update current position
    state.dx .= state.alpha.*state.s
    state.x .= state.x .+ state.dx
    retract!(method.manifold, state.x)

    lssuccess == false # break on linesearch error
end

function trace!(tr, d, state, iteration, method::BFGS, options, curr_time=time())
    dt = Dict()
    dt["time"] = curr_time
    if options.extended_trace
        dt["x"] = copy(state.x)
        dt["g(x)"] = copy(gradient(d))
        dt["~inv(H)"] = copy(state.invH)
        dt["Current step size"] = state.alpha
    end
    g_norm = norm(gradient(d), Inf)
    update!(tr,
    iteration,
    value(d),
    g_norm,
    dt,
    options.store_trace,
    options.show_trace,
    options.show_every,
    options.callback)
end

function update!(tr::OptimizationTrace{Tf},
              iteration::Integer,
              f_x::Tf,
              grnorm::Real,
              dt::Dict,
              store_trace::Bool,
              show_trace::Bool,
              show_every::Int = 1,
              callback = nothing) where {Tf}
    os = OptimizationState{Tf}(iteration, f_x, grnorm, dt)
    if store_trace
        push!(tr, os)
    end
    if show_trace
        if iteration % show_every == 0
            show(os)
            flush(stdout)
        end
    end
    if callback !== nothing && (iteration % show_every == 0)
        if store_trace
            stopped = callback(tr)
        else
            stopped = callback(os)
        end
    else
        stopped = false
    end
    stopped
end

function update_g!(d, state, method)
    # Update the function value and gradient
    value_gradient!(d, state.x)
    project_tangent!(method.manifold, gradient(d), state.x)
end

function update_h!(d, state, method::BFGS)
    n = length(state.x)
    # Measure the change in the gradient
    state.dg .= gradient(d) .- state.g_previous

    # Update the inverse Hessian approximation using Sherman-Morrison
    dx_dg = real(dot(state.dx, state.dg))
    if dx_dg > 0
        mul!(vec(state.u), state.invH, vec(state.dg))

        c1 = (dx_dg + real(dot(state.dg, state.u))) / (dx_dg' * dx_dg)
        c2 = 1 / dx_dg

        # invH = invH + c1 * (s * s') - c2 * (u * s' + s * u')
        if(state.invH isa Array) # i.e. not a CuArray
            invH = state.invH; dx = state.dx; u = state.u;
            @inbounds for j in 1:n
                c1dxj = c1 * dx[j]'
                c2dxj = c2 * dx[j]'
                c2uj  = c2 *  u[j]'
                for i in 1:n
                    invH[i, j] = muladd(dx[i], c1dxj, muladd(-u[i], c2dxj, muladd(c2uj, -dx[i], invH[i, j])))
                end
            end
        else
            mul!(state.invH,vec(state.dx),vec(state.dx)', c1,1)
            mul!(state.invH,vec(state.u ),vec(state.dx)',-c2,1)
            mul!(state.invH,vec(state.dx),vec(state.u )',-c2,1)
        end
    end
end

function optimize(f, g, initial_x::AbstractArray, method::AbstractOptimizer, options::Options)
    d = OnceDifferentiable1(f, g, initial_x, real(zero(eltype(initial_x))))
    state = initial_state(method, d, initial_x)

    t0 = time() # Initial time stamp used to control early stopping by options.time_limit
    tr = OptimizationTrace{typeof(value(d))}()
    tracing = options.store_trace || options.show_trace || options.extended_trace || options.callback !== nothing
    stopped, stopped_by_callback, stopped_by_time_limit = false, false, false
    f_limit_reached, g_limit_reached, h_limit_reached = false, false, false
    x_converged, f_converged, f_increased, counter_f_tol = false, false, false, 0

    f_converged, g_converged = initial_convergence(d, state, method, initial_x, options)
    converged = f_converged || g_converged
    # prepare iteration counter (used to make "initial state" trace entry)
    iteration = 0

    options.show_trace && print_header(method)
    _time = time()
    trace!(tr, d, state, iteration, method, options, _time-t0)
    ls_success::Bool = true
    while !converged && !stopped && iteration < options.iterations
        iteration += 1
        ls_success = !update_state!(d, state, method)
        if !ls_success
            break # it returns true if it's forced by something in update! to stop (eg dx_dg == 0.0 in BFGS, or linesearch errors)
        end
        update_g!(d, state, method) # TODO: Should this be `update_fg!`?
        x_converged, f_converged,
        g_converged, f_increased = assess_convergence(state, d, options)
        # For some problems it may be useful to require `f_converged` to be hit multiple times
        # TODO: Do the same for x_tol?
        counter_f_tol = f_converged ? counter_f_tol+1 : 0
        converged = x_converged || g_converged || (counter_f_tol > options.successive_f_tol)
        update_h!(d, state, method) # only relevant if not converged
        if tracing
            # update trace; callbacks can stop routine early by returning true
            stopped_by_callback = trace!(tr, d, state, iteration, method, options, time()-t0)
        end

        # Check time_limit; if none is provided it is NaN and the comparison
        # will always return false.
        _time = time()
        stopped_by_time_limit = _time-t0 > options.time_limit
        f_limit_reached = options.f_calls_limit > 0 && f_calls(d) >= options.f_calls_limit ? true : false
        g_limit_reached = options.g_calls_limit > 0 && g_calls(d) >= options.g_calls_limit ? true : false
        h_limit_reached = options.h_calls_limit > 0 && h_calls(d) >= options.h_calls_limit ? true : false

        if (f_increased && !options.allow_f_increases) || stopped_by_callback ||
            stopped_by_time_limit || f_limit_reached || g_limit_reached || h_limit_reached
            stopped = true
        end

        if g_calls(d) > 0 && !all(isfinite, gradient(d))
            @warn "Terminated early due to NaN in gradient."
            break
        end
        if h_calls(d) > 0 && !(d isa TwiceDifferentiableHV) && !all(isfinite, hessian(d))
            @warn "Terminated early due to NaN in Hessian."
            break
        end
    end # while

    # we can just check minimum, as we've earlier enforced same types/eltypes
    # in variables besides the option settings
    Tf = typeof(value(d))
    f_incr_pick = f_increased && !options.allow_f_increases
    stopped_by =(f_limit_reached=f_limit_reached,
                 g_limit_reached=g_limit_reached,
                 h_limit_reached=h_limit_reached,
                 time_limit=stopped_by_time_limit,
                 callback=stopped_by_callback,
                 f_increased=f_incr_pick)
    return MultivariateOptimizationResults(method,
                                        initial_x,
                                        pick_best_x(f_incr_pick, state),
                                        pick_best_f(f_incr_pick, state, d),
                                        iteration,
                                        iteration == options.iterations,
                                        x_converged,
                                        Tf(options.x_abstol),
                                        Tf(options.x_reltol),
                                        x_abschange(state),
                                        x_relchange(state),
                                        f_converged,
                                        Tf(options.f_abstol),
                                        Tf(options.f_reltol),
                                        f_abschange(d, state),
                                        f_relchange(d, state),
                                        g_converged,
                                        Tf(options.g_abstol),
                                        g_residual(d, state),
                                        f_increased,
                                        tr,
                                        f_calls(d),
                                        g_calls(d),
                                        h_calls(d),
                                        ls_success,
                                        options.time_limit,
                                        _time-t0,
                                        stopped_by,
                                        )
end

retract!(M::Manifold,x) = x
retract(M::Manifold,x) = retract!(M, copy(x))
project_tangent!(M::Manifold, g, x) = g

value(obj::AbstractObjective) = obj.F
function value(obj::ManifoldObjective)
    value(obj.inner_obj)
end

gradient(obj::AbstractObjective) = obj.DF
function gradient(obj::ManifoldObjective)
    gradient(obj.inner_obj)
end
function gradient!(obj::AbstractObjective, x)
    if x != obj.x_df
        gradient!!(obj, x)
    end
    gradient(obj)
end

function value_gradient!(obj::AbstractObjective, x)
    if x != obj.x_f && x != obj.x_df
        value_gradient!!(obj, x)
    elseif x != obj.x_f
        value!!(obj, x)
    elseif x != obj.x_df
        gradient!!(obj, x)
    end
    value(obj), gradient(obj)
end

function value_gradient!(obj::ManifoldObjective,x)
    xin = retract(obj.manifold, x)
    value_gradient!(obj.inner_obj,xin)
    project_tangent!(obj.manifold,gradient(obj.inner_obj),xin)
    return value(obj.inner_obj)
end

function value_gradient!!(obj::AbstractObjective, x)
    obj.f_calls .+= 1
    obj.df_calls .+= 1
    copyto!(obj.x_f, x)
    copyto!(obj.x_df, x)
    obj.F = obj.fdf(gradient(obj), x)
    value(obj), gradient(obj)
end

function _init_identity_matrix(x::AbstractArray{T}, scale::T = T(1)) where {T}
    x_ = reshape(x, :)
    Id = x_ .* x_' .* false
    idxs = diagind(Id)
    @. @view(Id[idxs]) = scale * true
    return Id
end

function perform_linesearch!(state, method, d)
    # Calculate search direction dphi0
    dphi_0 = real(dot(gradient(d), state.s))
    # reset the direction if it becomes corrupted
    if dphi_0 >= zero(dphi_0) && reset_search_direction!(state, d, method)
        dphi_0 = real(dot(gradient(d), state.s)) # update after direction reset
    end
    phi_0  = value(d)

    # Guess an alpha
    method.alphaguess!(method.linesearch!, state, phi_0, dphi_0, d)

    # Store current x and f(x) for next iteration
    state.f_x_previous = phi_0
    copyto!(state.x_previous, state.x)

    # Perform line search; catch LineSearchException to allow graceful exit
    try
        state.alpha, ϕalpha = method.linesearch!(d, state.x, state.s, state.alpha,
                               state.x_ls, phi_0, dphi_0)
        return true # lssuccess = true
    catch ex
        if isa(ex, LineSearchException)
            state.alpha = ex.alpha
            # We shouldn't warn here, we should just carry it to the output
            # @warn("Linesearch failed, using alpha = $(state.alpha) and
            # exiting optimization.\nThe linesearch exited with message:\n$(ex.message)")
            return false # lssuccess = false
        else
            rethrow(ex)
        end
    end
end

function print_header(method::AbstractOptimizer)
    @printf "Iter     Function value   Gradient norm \n"
end

function Base.show(io::IO, t::OptimizationState)
    @printf io "%6d   %14e   %14e\n" t.iteration t.value t.g_norm
    if !isempty(t.metadata)
        for (key, value) in t.metadata
            @printf io " * %s: %s\n" key value
        end
    end
    return
end

x_of_nans(x, Tf=eltype(x)) = fill!(Tf.(x), Tf(NaN))
alloc_DF(x, F::T) where T<:Number = x_of_nans(x, promote_type(eltype(x), T))

g_calls(r::MultivariateOptimizationResults) = r.g_calls
g_calls(d) = first(d.df_calls)
h_calls(r::MultivariateOptimizationResults) = r.h_calls
h_calls(d::OnceDifferentiable1) = 0
h_calls(d) = first(d.h_calls)
f_calls(d) = first(d.f_calls)

pick_best_x(f_increased, state) = f_increased ? state.x_previous : state.x
pick_best_f(f_increased, state, d) = f_increased ? state.f_x_previous : value(d)

function maxdiff(x::AbstractArray, y::AbstractArray)
    return mapreduce((a, b) -> abs(a - b), max, x, y)
end

f_abschange(d::AbstractObjective, state) = f_abschange(value(d), state.f_x_previous)
f_abschange(f_x::T, f_x_previous) where T = abs(f_x - f_x_previous)
f_abschange(r::MultivariateOptimizationResults) = r.f_abschange
f_relchange(d::AbstractObjective, state) = f_relchange(value(d), state.f_x_previous)
f_relchange(f_x::T, f_x_previous) where T = abs(f_x - f_x_previous)/abs(f_x)
f_relchange(r::MultivariateOptimizationResults) = r.f_relchange

x_abschange(r::MultivariateOptimizationResults) = r.x_abschange
x_relchange(r::MultivariateOptimizationResults) = r.x_relchange
x_abschange(state) = x_abschange(state.x, state.x_previous)
x_abschange(x, x_previous) = maxdiff(x, x_previous)
x_relchange(state) = x_relchange(state.x, state.x_previous)
x_relchange(x, x_previous) = maxdiff(x, x_previous)/maximum(abs, x)

g_residual(d, state) = g_residual(d)
g_residual(d::AbstractObjective) = g_residual(gradient(d))
g_residual(g) = maximum(abs, g)
g_residual(r::MultivariateOptimizationResults) = r.g_residual

function initial_convergence(d, state, method::AbstractOptimizer, initial_x, options)
    gradient!(d, initial_x)
    stopped = !isfinite(value(d)) || any(!isfinite, gradient(d))
    maximum(abs, gradient(d)) <= options.g_abstol, stopped
end

function converged(r::MultivariateOptimizationResults)
    conv_flags = r.x_converged || r.f_converged || r.g_converged
    x_isfinite = isfinite(x_abschange(r)) || isnan(x_relchange(r))
    f_isfinite = if r.iterations > 0
            isfinite(f_abschange(r)) || isnan(f_relchange(r))
        else
            true
        end
    g_isfinite = isfinite(g_residual(r))
    return conv_flags && all((x_isfinite, f_isfinite, g_isfinite))
end

# Default function for convergence assessment used by BFGSState
function assess_convergence(state::BFGSState, d, options::Options)
    x_converged, f_converged, f_increased, g_converged = false, false, false, false

    f_x = value(d)
    g_x = gradient(d)

    if x_abschange(state.x, state.x_previous) ≤ options.x_abstol
        x_converged = true
    end
    if x_abschange(state.x, state.x_previous) ≤ options.x_reltol * maximum(abs, state.x)
        x_converged = true
    end

    if f_abschange(f_x, state.f_x_previous) ≤ options.f_abstol
        f_converged = true
    end

    if f_abschange(f_x, state.f_x_previous) ≤ options.f_reltol*abs(f_x)
        f_converged = true
    end

    if f_x > state.f_x_previous
        f_increased = true
    end

    g_converged = g_residual(g_x) ≤ options.g_abstol

    return x_converged, f_converged, g_converged, f_increased
end
