#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/11/03
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

    res = optimize(𝑓, 𝐽!, 𝑎𝑏, max_iter = 500)
    
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

mutable struct BFGSState{Tx, Tm, T, G}
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

mutable struct BFGSOptimizationResults{Tx, Tc, Tf}
    initial_x::Tx
    minimizer::Tx
    minimum::Tf
    iterations::Int
    x_abschange::Tc
    x_relchange::Tc
    f_abschange::Tc
    f_relchange::Tc
    g_converged::Bool
    g_residual::Tc
end

include("hagerzhang.jl")

# Used for objectives and solvers where the gradient is available/exists
mutable struct BFGSDifferentiable{TF, TDF}
    ℱ! # objective, f
    𝒟! # (partial) derivative of objective, df
    𝐹 :: TF # cache for f output, F
    𝐷 :: TDF # cache for df output, DF
end

x_of_nans(x, Tf=eltype(x)) = fill!(Tf.(x), Tf(NaN))
alloc_DF(x, F::T) where T<:Number = x_of_nans(x, promote_type(eltype(x), T))

function BFGSDifferentiable(f, df, x::AbstractArray)
    F::Real = real(zero(eltype(x)))
    DF::AbstractArray = alloc_DF(x, F)
    BFGSDifferentiable(f, df, copy(F), copy(DF))
end

value(obj::BFGSDifferentiable) = obj.𝐹
gradient(obj::BFGSDifferentiable) = obj.𝐷
function value_gradient!(obj::BFGSDifferentiable, x)
    obj.𝒟!(gradient(obj), x)
    obj.𝐹 = obj.ℱ!(x)
end

function init_state(d::BFGSDifferentiable, initial_x::AbstractArray{T}) where T
    value_gradient!(d, initial_x)

    x_ = reshape(initial_x, :)
    invH0 = x_ .* x_' .* false
    idxs = diagind(invH0)
    scale = eltype(initial_x)(1)
    @. @view(invH0[idxs]) = scale * true

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

function update_state!(d::BFGSDifferentiable, state::BFGSState)
    T = eltype(state.s)
    # Set the search direction
    # Search direction is the negative gradient divided by the approximate Hessian
    mul!(vec(state.s), state.invH, vec(gradient(d)))
    rmul!(state.s, T(-1))

    # Maintain a record of the previous gradient
    copyto!(state.g_previous, gradient(d))

    # Determine the distance of movement along the search line
    # This call resets invH to initial_invH is the former in not positive
    # semi-definite
    lssuccess = perform_linesearch!(state, d)

    # Update current position
    state.dx .= state.alpha.*state.s
    state.x .= state.x .+ state.dx

    lssuccess == false # break on linesearch error
end

function trace!(d::BFGSDifferentiable, iteration, curr_time=time())
    dt = Dict()
    dt["time"] = curr_time
    g_norm = norm(gradient(d), Inf)

    if iteration % 1 == 0
        @printf("%6d   %14e   %14e\n", iteration, value(d), g_norm)
        if !isempty(dt)
            for (key, value) in dt
                @printf(" * %s: %s\n", key, value)
            end
        end
        flush(stdout)
    end
    false
end

# Update the function value and gradient
function update_g!(d::BFGSDifferentiable, state::BFGSState)
    value_gradient!(d, state.x)
end

function update_h!(d::BFGSDifferentiable, state::BFGSState)
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

function optimize(f, g, initial_x::AbstractArray; max_iter::I64 = 1000)
    d = BFGSDifferentiable(f, g, initial_x)
    state = init_state(d, initial_x)

    t0 = time() # Initial time stamp

    stopped = false

    g_converged = !isfinite(value(d)) || any(!isfinite, gradient(d)) #initial_convergence(d)
    # prepare iteration counter (used to make "initial state" trace entry)
    iteration = 0

    @printf "Iter     Function value   Gradient norm \n"

    _time = time()
    trace!(d, iteration, _time-t0)
    ls_success = true
    while !g_converged && !stopped && iteration < max_iter
        iteration += 1
        ls_success = !update_state!(d, state)
        if !ls_success
            break # it returns true if it's forced by something in update! to stop (eg dx_dg == 0.0 in BFGS, or linesearch errors)
        end
        update_g!(d, state)
        g_converged = (g_residual(d) ≤ 1e-8)
        update_h!(d, state) # only relevant if not converged

        # update trace
        trace!(d, iteration, time()-t0)

        _time = time()

        if !all(isfinite, gradient(d))
            @warn "Terminated early due to NaN in gradient."
            break
        end
    end # while

    # we can just check minimum, as we've earlier enforced same types/eltypes
    # in variables besides the option settings
    return BFGSOptimizationResults(initial_x,
                                        state.x,
                                        value(d),
                                        iteration,
                                        x_abschange(state),
                                        x_relchange(state),
                                        f_abschange(d, state),
                                        f_relchange(d, state),
                                        g_converged,
                                        g_residual(d)
    )
end

function perform_linesearch!(state::BFGSState, d::BFGSDifferentiable)
    # Calculate search direction dphi0
    dphi_0 = real(dot(gradient(d), state.s))
    # reset the direction if it becomes corrupted
    if dphi_0 >= zero(dphi_0)
        dphi_0 = real(dot(gradient(d), state.s)) # update after direction reset
    end
    phi_0  = value(d)

    # Guess an alpha
    guess = InitialStatic()
    linesearch = HagerZhang()
    guess(linesearch, state, phi_0, dphi_0, d)

    # Store current x and f(x) for next iteration
    state.f_x_previous = phi_0
    copyto!(state.x_previous, state.x)

    # Perform line search; catch LineSearchException to allow graceful exit
    try
        state.alpha, ϕalpha = linesearch(d, state.x, state.s, state.alpha,
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

function maxdiff(x::AbstractArray, y::AbstractArray)
    return mapreduce((a, b) -> abs(a - b), max, x, y)
end

f_abschange(d::BFGSDifferentiable, state) = abs(value(d) - state.f_x_previous)
f_relchange(d::BFGSDifferentiable, state) = abs(value(d) - state.f_x_previous)/abs(value(d))
x_abschange(state::BFGSState) = maxdiff(state.x, state.x_previous)
x_relchange(state::BFGSState) = maxdiff(state.x, state.x_previous)/maximum(abs, state.x)
g_residual(d::BFGSDifferentiable) = g_residual(gradient(d))
g_residual(g) = maximum(abs, g)

function converged(r::BFGSOptimizationResults)
    conv_flags = r.g_converged
    x_isfinite = isfinite(r.x_abschange) || isnan(r.x_relchange)
    f_isfinite = if r.iterations > 0
            isfinite(r.f_abschange) || isnan(r.f_relchange)
        else
            true
        end
    g_isfinite = isfinite(r.g_residual)
    return conv_flags && all((x_isfinite, f_isfinite, g_isfinite))
end
