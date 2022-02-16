#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/16
#

#=
*Remarks* : 

In ACFlow, the following kernrl functions are supported.

* Fermionic Imaginary-Time Kernel 

```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}} A(\omega)
\end{equation}
```

```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}}
\end{equation}
```

* Fermionic Matsubara Frequency Kernel

```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\epsilon
               \frac{1}{i\omega_n - \epsilon} A(\epsilon)
\end{equation}
```

```math
\begin{equation}
K(\omega_n,\epsilon) = \frac{1}{i\omega_n - \epsilon}
\end{equation}
```

* Bosonic Imaginary-Time Kernel

```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 - e^{-\beta\omega}} A(\omega)
\end{equation}
```

```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 - e^{-\beta\omega}}
\end{equation}
```

* Bosonic Matsubara Frequency Kernel

```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\epsilon
               \frac{1}{i\omega_n - \epsilon} A(\epsilon)
\end{equation}
```

```math
\begin{equation}
K(\omega_n,\epsilon) = \frac{1}{i\omega_n - \epsilon}
\end{equation}
```

* Bosonic Imaginary-Time Kernel (Symmetric Mode)

```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega} + e^{-(\beta - \tau)\omega}}
                      {2(1 - e^{-\beta\omega})}
\end{equation}
```

* Bosonic Matsubara Frequency Kernel (Symmetric Mode)

```math
\begin{equation}
K(\omega_n, \epsilon) = \frac{\epsilon^2}{\omega_n^2 + \epsilon^2}
\end{equation}
```
=#

"""
    build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)

Try to build kernel function in fermionic imaginary time axis.

See also: [`AbstractMesh`](@ref), [`FermionicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)
    ntime = fg.ntime
    nmesh = am.nmesh
    β = fg.β

    kernel = zeros(F64, ntime, nmesh)
    for i = 1:nmesh
        for j = 1:ntime
            kernel[j,i] = exp(-fg[j] * am[i]) / (1.0 + exp(-β * am[i]))
        end
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)

Try to build kernel function in fermionic Matsubara frequency axis.

See also: [`AbstractMesh`](@ref), [`FermionicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = fg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    if blur > 0.0
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)
        Mx = reshape(gaussian, (1, 1, nsize))
        Mg = reshape(fg.ω, (nfreq, 1, 1))
        Mm = reshape(am.mesh, (1, nmesh, 1))
        Mb = reshape(bmesh, (1, 1, nsize))

        integrand = Mx ./ (im * Mg .- Mm .- Mb)

        for i = 1:nmesh
            for j = 1:nfreq
                _kernel[j,i] = simpson(bmesh, integrand[j,i,:])
            end
        end
    else
        for i = 1:nmesh
            for j = 1:nfreq
                _kernel[j,i] = 1.0 / (im * fg[j] - am[i])
            end
        end
    end

    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)

Try to build kernel function in bosonic imaginary time axis.

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = fg.β

    kernel = zeros(F64, ntime, nmesh)
    for i = 1:nmesh
        for j = 1:ntime
            kernel[j,i] = exp(-bg[j] * am[i]) / (1.0 - exp(-β * am[i]))
        end
    end
    @. kernel[:,1] = 1.0

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build kernel function in bosonic Matsubara frequency axis.

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    nfreq = bg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    for i = 1:nmesh
        for j = 1:nfreq
            _kernel[j,i] = 1.0 / (im * bg[j] - am[i])
        end
    end

    if am[1] == 0.0 && bg[1] == 0.0
        _kernel[1,1] = 1.0
    end

    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)

Try to build kernel function in bosonic imaginary time axis.

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    kernel = zeros(F64, ntime, nmesh)
    for i = 1:nmesh
        r = 0.5 / (1.0 - exp(-β * am[i]))
        for j = 1:ntime
            kernel[j,i] = r * (exp(-am[i] * bg[j]) + exp(-am[i] * (β - bg[j])))
        end
    end
    @. kernel[:,1] = 1.0

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build kernel function in bosonic Matsubara frequency axis.

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = bg.nfreq
    nmesh = am.nmesh

    kernel = zeros(F64, nfreq, nmesh)

    if blur > 0.0
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)
        Mx = reshape(gaussian, (1, 1, nsize))
        Mg = reshape(bg.ω, (nfreq, 1, 1))
        Mm = reshape(am.mesh, (1, nmesh, 1))
        Mb = reshape(bmesh, (1, 1, nsize))

        integrand_1 = Mx .* ((Mb .+ Mm) .^ 2.0) ./ ((Mb .+ Mm) .^ 2.0 .+ Mg .^ 2.0)
        integrand_2 = Mx .* ((Mb .- Mm) .^ 2.0) ./ ((Mb .- Mm) .^ 2.0 .+ Mg .^ 2.0)
        for j = 1:nmesh
            integrand_1[1,j,:] .= gaussian
            integrand_2[1,j,:] .= gaussian
        end
        integrand = (integrand_1 + integrand_2) / 2.0
        for i = 1:nmesh
            for j = 1:nfreq
                kernel[j,i] = simpson(bmesh, integrand[j,i,:])
            end
        end
    else
        for i = 1:nmesh
            for j = 1:nfreq
                kernel[j,i] = am[i] ^ 2.0 / ( bg[j] ^ 2.0 + am[i] ^ 2.0 )
            end
        end
        if am[1] == 0.0 && bg[1] == 0.0
            kernel[1,1] = 1.0
        end
    end

    return kernel
end

"""
    make_blur(am::AbstractMesh, A::Vector{F64}, blur::F64)

Try to blur the given spectrum `A`, which is defined in `am`. And `blur`
is the blur parameter.
"""
function make_blur(am::AbstractMesh, A::Vector{F64}, blur::F64)
    ktype = get_c("ktype")

    spl = nothing
    if ktype == "fermi" || ktype == "boson"
        spl = CubicSpline(A, am.mesh)
    else
        vM = vcat(-am.mesh[end:-1:2], am.mesh)
        vA = vcat(A[end:-1:2], A)
        spl = CubicSpline(vA, vM)
    end

    bmesh, gaussian = make_gauss_peaks(blur)

    nsize = length(bmesh)
    nmesh = am.nmesh

    Mb = reshape(bmesh, (1, nsize))
    Mx = reshape(gaussian, (1, nsize))
    Mm = reshape(am.mesh, (nmesh, 1))
    integrand = Mx .* spl.(Mm .+ Mb)

    for i = 1:nmesh
        A[i] = simpson(bmesh, integrand[i,:])
    end
end

"""
    make_singular_space(kernel::Matrix{F64})

Perform singular value decomposition for the input matrix.
"""
function make_singular_space(kernel::Matrix{F64})
    U, S, V = svd(kernel)
    n_svd = count(x -> x ≥ 1e-10, S)
    U_svd = U[:,1:n_svd]
    V_svd = V[:,1:n_svd]
    S_svd = S[1:n_svd]

    return U_svd, V_svd, S_svd
end

"""
    make_gauss_peaks(blur::F64)

Try to generate a series of gaussian peaks along a linear mesh, whose
energy range is `[-5 * blur, +5 * blur]`.
"""
function make_gauss_peaks(blur::F64)
    @assert blur > 0.0
    nsize = 201
    bmesh = collect(LinRange(-5.0 * blur, 5.0 * blur, nsize))
    norm = 1.0 / (blur * sqrt(2.0 * π))
    gaussian = norm * exp.(-0.5 * (bmesh / blur) .^ 2.0)
    return bmesh, gaussian
end
