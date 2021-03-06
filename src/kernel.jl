#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/05/05
#

#=
*Remarks* :

The purpose of `ACFlow` is to solve the following equation:

```math
\begin{equation}
\mathbf{G} = \mathbf{KA}.
\end{equation}
```

Here, ``\mathbf{G}``, ``\mathbf{K}``, and ``\mathbf{A}`` are input
green's function, kernel function, and spectral density, respectively.
`ACFlow` supports various kernel functions. They are summaried as follows.

**For Fermionic Green's Functions**

In imaginary-time axis, we have

```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}} A(\omega),
\end{equation}
```

```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}}.
\end{equation}
```

In Matsubara frequency axis, we have

```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{1}{i\omega_n - \omega} A(\omega),
\end{equation}
```

```math
\begin{equation}
K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega},
\end{equation}
```

where ``\omega_n`` is a Matsubara frequencies equal to ``(2n + 1)\pi/\beta``.

These kernel functions are for the finite temperature Green's function of
fermions,

```math
G(\tau) = -\langle \mathcal{T} c(\tau) c^{\dagger}(0)\rangle.
```

``G(\tau)`` must fulfil the anti-periodicity condition,

```math
G(\tau + \beta) = -G(\tau).
```

And its Fourier transform and inverse Fourier transform are given by

```math
G(i\omega_n) = \int^{\beta}_0 d\tau\ e^{-i\omega_n \tau} G(\tau),
```

```math
G(\tau) = \frac{1}{\beta} \sum_n e^{i\omega_n \tau} G(i\omega_n).
```

It is possible to analytically continue similar anti-periodic
functions, such as fermionic self-energy function ``\Sigma``,
with these kernel functions. For the self-energies, it is addtionally
required that the constant contribution ``\Sigma(i\infty)`` is
subtracted from ``\Sigma(i\omega_n)``.

**For Bosonic Green's Functions**

In imaginary-time axis, we have

```math
\begin{equation}
G_{B}(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}
          \frac{A(\omega)}{\omega},
\end{equation}
```

```math
\begin{equation}
K(\tau,\omega) = \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}.
\end{equation}
```

if ``\omega = 0``,

```math
\begin{equation}
K(\tau,0) = \frac{1}{\beta}.
\end{equation}
```

In Matsubara frequency axis, we have

```math
\begin{equation}
G_{B}(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{\omega}{i\omega_n - \omega}
               \frac{A(\omega)}{\omega},
\end{equation}
```

```math
\begin{equation}
K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega},
\end{equation}
```

where ``\omega_n`` is a Matsubara frequencies equal to ``2n\pi/\beta``.
``G_{B}(\tau)`` is related to ``G_{B}(i\omega_n)`` via the Fourier
transform as mentioned above.

Especially, if ``\omega_n = 0`` and ``\omega = 0``,

```math
\begin{equation}
K(0,0) = -1.
\end{equation}
```

These kernel functions are for finite temperature correlation function of
boson-like operators ``B`` and ``B^{\dagger}``,

```math
\chi_{B}(\tau) = \langle \mathcal{T} B(\tau) B^{\dagger}(0)\rangle.
```

``\chi_{B}(\tau)`` must be ``\beta``-periodic, i.e.,

```math
\chi_{B}(\tau + \beta) = \chi_{B}(\tau).
```

Typical examples of such functions are Green's function of bosons

```math
G_{b}(\tau) = \langle \mathcal{T} b(\tau) b^{\dagger}(0)\rangle,
```

and the transverse spin susceptibility

```math
\chi_{+-}(\tau) = \langle \mathcal{T} S_{+}(\tau) S_{-}(0) \rangle.
```

**For Correlator of Hermitian Operator**

In imaginary-time axis, we have

```math
\begin{equation}
K(\tau,\omega) = \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]}
                      {2(1 - e^{-\beta\omega})}.
\end{equation}
```

if ``\omega = 0``,

```math
K(\tau,0) = \frac{1}{\beta}.
```

In Matsubara frequency axis, we have

```math
\begin{equation}
K(\omega_n, \omega) = \frac{\omega^2}{\omega_n^2 + \omega^2}.
\end{equation}
```

Especially, if ``\omega_n = 0`` and ``\omega = 0``,

```math
\begin{equation}
K(0,0) = 1.
\end{equation}
```

This is a special case of the previous observable kind with ``B = B^{\dagger}``
, and its use is in general prefered due to the reduced ``A(\omega)``
definition domain. The most widely used observables of this kind are
the longitudinal spin susceptibility,

```math
\chi_{zz}(\tau) = \langle S_z(\tau) S_z(0) \rangle,
```

and the charge susceptibility,

```math
\chi(\tau) = \langle N(\tau) N(0) \rangle.
```

**Cautions:**

For bosonic-like system, the calculated spectral density might be not
``A(\omega)``. It could be ``A(\omega)/\omega`` or any others. It in
fact depends on the type of the correlation function and how the kernel
function is defined.
=#

"""
    build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)

Try to build fermionic kernel function in imaginary time axis.

See also: [`AbstractMesh`](@ref), [`FermionicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)
    ntime = fg.ntime
    nmesh = am.nmesh
    ?? = fg.??

    kernel = zeros(F64, ntime, nmesh)
    for i = 1:nmesh
        for j = 1:ntime
            kernel[j,i] = exp(-fg[j] * am[i]) / (1.0 + exp(-?? * am[i]))
        end
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)

Try to build fermionic kernel function in Matsubara frequency axis. This
function support preblur algorithm.

See also: [`AbstractMesh`](@ref), [`FermionicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = fg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                _kernel[j,i] = 1.0 / (im * fg[j] - am[i])
            end
        end
    else
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)
        integrand = zeros(C64, nsize)
        for i = 1:nmesh
            for j = 1:nfreq
                z = im * fg[j] - am[i]
                for k = 1:nsize
                    integrand[k] = gaussian[k] / (z - bmesh[k])
                end
                _kernel[j,i] = simpson(bmesh, integrand)
            end
        end
    end

    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)

Try to build bosonic kernel function in imaginary time axis.

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    ?? = bg.??

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        for j = 1:ntime
            kernel[j,i] = am[i] * exp(-bg[j] * am[i]) / (1.0 - exp(-?? * am[i]))
        end
    end
    #
    @. kernel[:,1] = 1.0 / ??

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis.

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    nfreq = bg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    for i = 1:nmesh
        for j = 1:nfreq
            _kernel[j,i] = am[i] / (im * bg[j] - am[i])
        end
    end

    if am[1] == 0.0 && bg[1] == 0.0
        _kernel[1,1] = -1.0
    end

    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)

Try to build bosonic kernel function in imaginary time axis (just for
correlator of Hermitian operator only).

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    ?? = bg.??

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        r = am[i] * 0.5 / (1.0 - exp(-?? * am[i]))
        for j = 1:ntime
            kernel[j,i] = r * (exp(-am[i] * bg[j]) + exp(-am[i] * (?? - bg[j])))
        end
    end
    #
    @. kernel[:,1] = 1.0 / ??

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis (just
for correlator of Hermitian operator only). This function support preblur
algorithm.

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = bg.nfreq
    nmesh = am.nmesh

    kernel = zeros(F64, nfreq, nmesh)

    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                kernel[j,i] = am[i] ^ 2.0 / ( bg[j] ^ 2.0 + am[i] ^ 2.0 )
            end
        end
        if am[1] == 0.0 && bg[1] == 0.0
            kernel[1,1] = 1.0
        end
    else
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)

        I??? = zeros(F64, nsize)
        I??? = zeros(F64, nsize)
        I??? = zeros(F64, nsize)

        for i = 1:nmesh
            for j = 1:nfreq
                g?? = bg[j] ^ 2.0
                for k = 1:nsize
                    A?? = (bmesh[k] + am[i]) ^ 2.0; I???[k] = A?? / (A?? + g??)
                    B?? = (bmesh[k] - am[i]) ^ 2.0; I???[k] = B?? / (B?? + g??)
                end
                if i == 1 && j == 1
                    I??? .= 1.0
                    I??? .= 1.0
                end
                @. I??? = (I??? + I???) * gaussian / 2.0
                kernel[j,i] = simpson(bmesh, I???)
            end
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
    ktype = get_b("ktype")

    spl = nothing
    if ktype == "fermi" || ktype == "boson"
        spl = CubicSplineInterpolation(A, am.mesh)
    else
        vM = vcat(-am.mesh[end:-1:2], am.mesh)
        vA = vcat(A[end:-1:2], A)
        spl = CubicSplineInterpolation(vA, vM)
    end

    bmesh, gaussian = make_gauss_peaks(blur)

    nsize = length(bmesh)
    nmesh = length(am)

    Mb = reshape(bmesh, (nsize, 1))
    Mx = reshape(gaussian, (nsize, 1))
    Mm = reshape(am.mesh, (1, nmesh))
    I = Mx .* spl.(Mm .+ Mb)

    for j = 1:nmesh
        A[j] = simpson(bmesh, view(I, :, j))
    end
end

"""
    make_singular_space(kernel::Matrix{F64})

Perform singular value decomposition for the input matrix.
"""
function make_singular_space(kernel::Matrix{F64})
    U, S, V = svd(kernel)

    n_svd = count(x -> x ??? 1e-10, S)
    U_svd = U[:,1:n_svd]
    V_svd = V[:,1:n_svd]
    S_svd = S[1:n_svd]

    return U_svd, V_svd, S_svd
end

#=
*Remarks* : About Preblur

The Gaussian is

```math
\begin{equation}
g(x) = \frac{1}{b \sqrt{2 \pi}}
       \exp\left(-\frac{x^2}{2b^2}\right).
\end{equation}
```

Here `b` is the blur parameter. In the fermionic case, the convolution
can be written as:

```math
\begin{equation}
K_{preblur}(i\nu_n,\omega) = \int_{-5b}^{5b} dx~
    \frac{g(x)}{i\nu_n - x - \omega}.
\end{equation}
```

In the bosonic case, the convolution can be written as

```math
K_{preblur}(i\omega_n,\nu) = \frac{1}{2} \int_{-5b}^{5b} dx~g(x)
    \left[
        \frac{(x+\nu)^2 }{(x+\nu)^2 + \omega_n^2} +
        \frac{(x-\nu)^2 }{(x-\nu)^2 + \omega_n^2}
    \right].
```

Integration over the Gaussian from \(-5b\) to \(5b\) is certainly sufficient.
=#

"""
    make_gauss_peaks(blur::F64)

Try to generate a series of gaussian peaks along a linear mesh, whose
energy range is `[-5 * blur, +5 * blur]`.
"""
function make_gauss_peaks(blur::F64)
    @assert blur > 0.0
    nsize = 201
    bmesh = collect(LinRange(-5.0 * blur, 5.0 * blur, nsize))
    norm = 1.0 / (blur * sqrt(2.0 * ??))
    gaussian = norm * exp.(-0.5 * (bmesh / blur) .^ 2.0)
    return bmesh, gaussian
end
