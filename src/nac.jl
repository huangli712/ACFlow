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
        #J .= gradient(𝑓, x)[1]
        J .= finite_difference_gradient(𝑓, x)
        @show J
 #       exit()
    end

    function fun(x::Vector{C64}, nac::NevanACContext, ℋ::Array{APC,2})
        f = smooth_norm(nac, ℋ, x)
        J = gradient(𝑓, x)[1]
        return f, J
    end

    #solution, call = newton(fun, 𝑎𝑏, nac, ℋ)
    #@show solution, call

    res = optimize(𝑓, 𝐽!, 𝑎𝑏, BFGS(), 
                   Optim.Options(iterations = 500, show_trace = true))
    #@show Optim.minimizer(res)

    #exit()
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(ℋ, Optim.minimizer(res))

    if causality && (Optim.converged(res))
        nac.hopt = H
        nac.𝑎𝑏 = Optim.minimizer(res)
        nac.ℋ = ℋ
    end
    
    return causality, (Optim.converged(res))
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

@inline function compute_epsilon(x::T, relstep::Real, absstep::Real) where {T<:Number}
    return max(relstep*abs(x), absstep)
end

@inline function default_relstep(::Type{T}) where {T<:Number}
    return cbrt(eps(real(T)))
end

struct GradientCache{CacheType1,CacheType2,CacheType3,CacheType4}
    c1::CacheType2
    c3::CacheType4
end

function GradientCache(x)
    typeof(x) <: AbstractArray # the vector->scalar case
    _c1 = nothing
    _c2 = nothing
    _c3 = zero(x)

    GradientCache{Nothing,typeof(_c1),typeof(_c2),typeof(_c3)}(_c1, _c3)
end

function finite_difference_gradient(f, x::StridedVector{<:Number})
    df = zero(eltype(x)) .* x
    cache = GradientCache(x)
    relstep=default_relstep(eltype(x))
    absstep=relstep

    # c1 is x1 if we need a complex copy of x, otherwise Nothing
    # c2 is Nothing
    c1, c3 = cache.c1, cache.c3

    fdtype=Val(:central)
    if fdtype != Val(:complex)
        if eltype(df) <: Complex && !(eltype(x) <: Complex)
            copyto!(c1, x)
        end
    end
    copyto!(c3, x)
    
    @inbounds for i ∈ eachindex(x)
        epsilon = compute_epsilon(x[i], relstep, absstep)
        x_old = x[i]
        c3[i] += epsilon
        dfi = f(c3)
        c3[i] = x_old - epsilon
        dfi -= f(c3)
        c3[i] = x_old
        df[i] = real(dfi / (2 * epsilon))
        if eltype(df) <: Complex
            if eltype(x) <: Complex
                c3[i] += im * epsilon
                dfi = f(c3)
                c3[i] = x_old - im * epsilon
                dfi -= f(c3)
                c3[i] = x_old
            else
                c1[i] += im * epsilon
                dfi = f(c1)
                c1[i] = x_old - im * epsilon
                dfi -= f(c1)
                c1[i] = x_old
            end
            df[i] -= im * imag(dfi / (2 * im * epsilon))
        end
    end
    df
end
