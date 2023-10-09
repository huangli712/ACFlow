#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/10/09
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
* Hopt -> Optimal value of H, the order of Hardy basis functions.
"""
mutable struct NevanACContext
    Gᵥ   :: Vector{APC}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    Φ    :: Vector{APC}
    𝒜    :: Array{APC,3}
    ℋ    :: Array{APC,2}
    𝑎𝑏   :: Vector{C64}
    Hopt :: I64
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

    return NevanACContext(Gᵥ, grid, mesh, Φ, 𝒜, ℋ, 𝑎𝑏, 1)
end

"""
    run(nac::NevanACContext)

Perform Hardy basis optimization to smooth the spectrum.
"""
function run(nac::NevanACContext)
    hardy = get_n("hardy")
    if hardy
        calc_hoptim(nac)
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
    # Note that _G is actually 𝑁G, so there is a `-` symbol.
    _G = C64.(calc_green(nac.𝒜, nac.ℋ, nac.𝑎𝑏))
    fwrite && write_complete(nac.mesh, -_G)

    # Calculate and write the spectral function
    Aout = F64.(imag.(_G) ./ π)
    fwrite && write_spectrum(nac.mesh, Aout)

    # Regenerate the input data and write them
    kernel = make_kernel(nac.mesh, nac.grid)
    G = reprod(nac.mesh, kernel, Aout)
    fwrite && write_backward(nac.grid, G)

    return Aout, _G
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

=#

"""
    precompute(grid::AbstractGrid,
               mesh::AbstractMesh,
               Gᵥ::Vector{APC})

Precompute some key quantities. Here `Gᵥ` is input data, `grid` is the
grid for input data, and `mesh` is the mesh for output spectrum.
"""
function precompute(grid::AbstractGrid,
                    mesh::AbstractMesh,
                    Gᵥ::Vector{APC})
    # Evaluate ϕ and `abcd` matrices
    Φ = calc_phis(grid, Gᵥ)
    𝒜 = calc_abcd(grid, mesh, Φ)

    # Allocate memory for evaluating θ
    ℋ = calc_hmatrix(mesh, 1)
    𝑎𝑏 = zeros(C64, 2)

    return Φ, 𝒜, ℋ, 𝑎𝑏
end

"""
    calc_mobius(z::Vector{APC})

A direct Mobius transformation.
"""
function calc_mobius(z::Vector{APC})
    _z = similar(z)
    @. _z = (z - im) / (z + im)
    return _z
end

"""
    calc_inv_mobius(z::Vector{APC})

An inverse Mobius transformation.
"""
function calc_inv_mobius(z::Vector{APC})
    _z = similar(z)
    @. _z = im * (one(APC) + z) / (one(APC) - z)
    return _z
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

    # Then we evaluate θ according Eq.(7)
    num = 𝒜[1,1,:] .* θₘ₊₁ .+ 𝒜[1,2,:]
    den = 𝒜[2,1,:] .* θₘ₊₁ .+ 𝒜[2,2,:]
    θ = num ./ den

    return θ
end

"""
    calc_green(𝒜::Array{APC,3}, ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})

θ is back transformed to a Nevanlinna interpolant via the inverse Mobius
transform. Here, `𝒜` (`abcd` matrix), `ℋ` (Hardy matrix), and `𝑎𝑏` are
used to evaluate θ.
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
    calc_hoptim(sol::NevanACContext)

Try to perform Hardy basis optimization. Such that the Hardy matrix ℋ
and the corresponding coefficients 𝑎𝑏 are updated. They are used to
calculate θ, which is then back transformed to generate smooth G (i.e.,
the spectrum) at real axis.
"""
function calc_hoptim(sol::NevanACContext)
    hmax = get_n("hmax")

    h = 1
    while h ≤ hmax
        println("H = $h")

        causality, optim = hardy_optim!(sol, h)

        # break if we find optimal H in which causality is preserved
        # and optimize is successful
        if causality && optim
            break
        else
            h = h + 1
        end
    end
end

function calc_loss(sol::NevanACContext, 𝑎𝑏::Vector{C64}, ℋ::Array{APC,2})
    #theta = calc_theta(sol.𝒜, ℋ, 𝑎𝑏)

    #green = im * (one(APC) .+ theta) ./ (one(APC) .- theta)
    green = calc_green(sol.𝒜, ℋ, 𝑎𝑏)
    A = F64.(imag(green)./pi)

    tot_int = trapz(sol.mesh, A)
    second_der = integrate_squared_second_deriv(sol.mesh.mesh, A) 

    alpha = get_n("alpha")
    func = abs(1.0-tot_int)^2 + alpha*second_der

    return func
end

function hardy_optim!(sol::NevanACContext, H::I64)::Tuple{Bool, Bool}
    ℋₗ = calc_hmatrix(sol.mesh, H)
    𝑎𝑏 = zeros(C64, 2*H)

    function functional(x::Vector{C64})::F64
        return calc_loss(sol, x, ℋₗ)
    end

    function jacobian(J::Vector{C64}, x::Vector{C64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, 𝑎𝑏, BFGS(), 
                   Optim.Options(iterations = 500, show_trace = true))
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(ℋₗ, Optim.minimizer(res))

    if causality && (Optim.converged(res))
        sol.Hopt = H
        sol.𝑎𝑏 = Optim.minimizer(res)
        sol.ℋ = ℋₗ
    end
    
    return causality, (Optim.converged(res))
end

"""
Compute second derivative
If the length of `x` and `y` is `N`, the length of the returned vector is `N-2`.
"""
function second_deriv(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x and y must be the same length")) 
    end

    N = length(x)
    dx_backward = view(x, 2:(N-1)) - view(x, 1:(N-2))
    dx_forward = view(x, 3:N) - view(x, 2:(N-1))

    y_forward = view(y, 3:N)
    y_mid = view(y, 2:(N-1))
    y_backward = view(y, 1:(N-2))

    n = dx_backward .* y_forward + dx_forward .* y_backward - (dx_forward + dx_backward) .* y_mid
    d = (dx_forward.^2) .* dx_backward + (dx_backward.^2) .* dx_forward
    return 2 .* n ./ d
end

"""
Integrate the squarre of the abs of the second derivative
"""
function integrate_squared_second_deriv(x::AbstractVector, y::AbstractVector)
    N = length(x)
    sd = second_deriv(x, y)

    x_sd = view(x, 2:(N-1))
    return trapz(x_sd, abs.(sd) .^ 2)
end

function test_pick(wn::Vector{APC}, gw::Vector{APC}, Nopt::I64)
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

function check_causality(ℋ::Array{APC,2}, 𝑎𝑏::Vector{C64})
    param = ℋ * 𝑎𝑏
    @show typeof(param), size(ℋ), size(𝑎𝑏), size(param)

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        println("max_theta=",max_theta)
        println("hardy optimization was success.")
        causality = true
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end

    return causality
end