#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/14
#

#=
*Remarks* :

**Preliminary Knowledge: Background**

The purpose of `ACFlow` is to solve the following equation:

```math
\begin{equation}
\mathbf{G} = \mathbf{KA}.
\end{equation}
```

Here, ``\mathbf{G}``, ``\mathbf{K}``, and ``\mathbf{A}`` are the input
Green's function, kernel function, and spectral density, respectively.
`ACFlow` supports various kernel functions. Here, we would like to dive
into them.

---

**Preliminary Knowledge: Imaginary-Time Green's Function**

The single-particle imaginary time Green's function reads,

```math
\begin{equation}
G_{F/B}(\tau) = \langle \mathcal{T}_{\tau} c(\tau) c^{\dagger}(0)\rangle,
\end{equation}
```

where `F` denotes fermions and `B` denotes bosons.

For fermions, ``G_{F}(\tau)`` must fulfil the anti-periodicity condition,

```math
\begin{equation}
G_{F}(\tau + \beta) = -G_{F}(\tau).
\end{equation}
```

For bosons, ``G_{B}(\tau)`` must be ``\beta``-periodic, i.e.,

```math
\begin{equation}
G_{B}(\tau + \beta) = G_{B}(\tau).
\end{equation}
```

---

**Preliminary Knowledge: Fourier Transformation**

The imaginary time Green's function ``G(\tau)`` and Matsubara Green's
function ``G(i\omega_n)`` are connected by the Fourier transformation
and inverse Fourier transformation:

```math
\begin{equation}
G(i\omega_n) = \int^{\beta}_0 d\tau\ e^{-i\omega_n \tau} G(\tau),
\end{equation}
```

```math
\begin{equation}
G(\tau) = \frac{1}{\beta} \sum_n e^{i\omega_n \tau} G(i\omega_n),
\end{equation}
```

where ``\beta`` is the inverse temperature (``\beta = 1/T``), ``\tau``
denotes the imaginary time and ``\omega_n`` means the Matsubara frequency.
Note that ``\omega_n = (2n + 1) \pi / \beta`` and ``2n \pi / \beta`` for
fermions and bosons, respectively.

---

**For Fermionic Green's Functions**

In imaginary time axis, we have

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

The spectral density ``A(\omega)`` is defined on ``(-\infty,\infty)``.
It is causal, i.e., ``A(\omega) \ge 0``.

It is also possible to analytically continue similar anti-periodic
functions, such as the fermionic self-energy function ``\Sigma(i\omega_n)``,
with these kernel functions. For the self-energies, it is addtionally
required that the constant contribution ``\Sigma(i\infty)`` is
subtracted from ``\Sigma(i\omega_n)`` beforehand.

---

**For Bosonic Green's Functions**

In imaginary time axis, we have

```math
\begin{equation}
G_{B}(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 - e^{-\beta\omega}}
          A(\omega),
\end{equation}
```

For bosonic system, the spectral density obey the following constraint:

```math
\begin{equation}
\text{sign}(\omega) A(\omega) \ge 0.
\end{equation}
```

It is quite convenient to introduce a new variable ``\tilde{A}(\omega)``:

```math
\begin{equation}
\tilde{A}(\omega) = \frac{A(\omega)}{\omega}.
\end{equation}
```

Clearly, ``\tilde{A}(\omega) \ge 0``. Thus, we have

```math
\begin{equation}
G_{B}(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}
          \tilde{A}(\omega),
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
               \frac{1}{i\omega_n - \omega} A(\omega)
                 = \int^{+\infty}_{-\infty} d\omega
               \frac{\omega}{i\omega_n - \omega} \tilde{A}(\omega)
\end{equation}
```

```math
\begin{equation}
K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega},
\end{equation}
```

Especially, if ``\omega_n = 0`` and ``\omega = 0``,

```math
\begin{equation}
K(0,0) = -1.
\end{equation}
```

Here, the spectral density ``A(\omega)``, or equivalently
``\tilde{A}(\omega)``, is defined on ``(-\infty,\infty)``.

Typical examples of this case include Green's function of bosons

```math
\begin{equation}
G_{b}(\tau) = \langle \mathcal{T}_{\tau} b(\tau) b^{\dagger}(0)\rangle,
\end{equation}
```

and the transverse spin susceptibility

```math
\begin{equation}
\chi_{+-}(\tau) = \langle \mathcal{T}_{\tau} S_{+}(\tau) S_{-}(0) \rangle.
\end{equation}
```

---

**For Correlator of Hermitian Operator**

In imaginary time axis, we have

```math
\begin{equation}
G_{B}(\tau) = \int^{\infty}_{0} d\omega
              \frac{e^{-\tau\omega} + e^{-(\beta - \tau)\omega}}
                   {1 - e^{-\beta\omega}}
              A(\omega).
\end{equation}
```

We introduce ``\tilde{A}(\omega) = A(\omega)/\omega`` again:

```math
\begin{equation}
G_{B}(\tau) = \int^{\infty}_{0} d\omega
              \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]}
                   {1 - e^{-\beta\omega}}
              \tilde{A}(\omega).
\end{equation}
```

```math
\begin{equation}
K(\tau,\omega) =
    \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]}
    {1 - e^{-\beta\omega}}.
\end{equation}
```

if ``\omega = 0``,

```math
\begin{equation}
K(\tau,0) = \frac{2}{\beta}.
\end{equation}
```

In Matsubara frequency axis, we have

```math
\begin{equation}
G_{B}(i\omega_n) = \int^{\infty}_{0} d\omega
                   \frac{-2\omega^2}{\omega_n^2 + \omega^2} \tilde{A}(\omega).
\end{equation}
```

```math
\begin{equation}
K(\omega_n, \omega) = \frac{-2\omega^2}{\omega_n^2 + \omega^2}.
\end{equation}
```

Especially, if ``\omega_n = 0`` and ``\omega = 0``,

```math
\begin{equation}
K(0,0) = -2.
\end{equation}
```

This is a special case of the previous observable kind with ``b = b^{\dagger}``
, and its use is in general preferred due to the reduced ``A(\omega)``
definition domain. Here, the spectral density ``A(\omega)``, or equivalently
``\tilde{A}(\omega)``, is defined on ``(0,\infty)``. Note that ``A(\omega)``
is an odd function, hence ``\tilde{A}(\omega)`` is an even function. The most
widely used observables of this kind are the longitudinal spin susceptibility,

```math
\begin{equation}
\chi_{zz}(\tau) = \langle S_z(\tau) S_z(0) \rangle,
\end{equation}
```

and the charge susceptibility,

```math
\begin{equation}
\chi(\tau) = \langle N(\tau) N(0) \rangle.
\end{equation}
```

**Cautions:**

For bosonic-like system, the calculated spectral density might be not
``A(\omega)``. It could be ``A(\omega)/\omega`` or any others. It in
fact depends on the type of the correlation function and how the kernel
function is defined.

**Reference:**

[1] J. Schott, *et al.*, Phys. Rev. B **94**, 245140 (2016).

[2] O. Gunnarsson, *et al.*, Phys. Rev. B **82**, 165125 (2010).

[3] M. Jarrell, *et al.*, Phys. Rep. **269**, 133 (1996).

---

=#

"""
    build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)

Try to build fermionic kernel function in imaginary time axis.

### Arguments
* am -> Real frequency mesh.
* fg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`FermionicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)
    ntime = fg.ntime
    nmesh = am.nmesh
    β = fg.β

    kernel = zeros(F64, ntime, nmesh)

    #
    # Old implementation
    # exp(βω ≈ 700.0) will throw an Inf and K becomes NaN.
    # We should avoid this situation.
    #
    # for i = 1:nmesh
    #     de = 1.0 + exp(-β * am[i])
    #     for j = 1:ntime
    #         kernel[j,i] = exp(-fg[j] * am[i]) / de
    #     end
    # end
    #
    # New implementation to avoid numerical overflow
    #
    for i = 1:nmesh
        # ω ≥ 0.0
        if am[i] ≥ 0.0
            de = 1.0 + exp(-β * am[i])
            for j = 1:ntime
                kernel[j,i] = exp(-fg[j] * am[i]) / de
            end
        # ω < 0.0
        else
            for j = 1:ntime
                # Now denominator goes to infinity and K goes to zero.
                if (fg[j] - β) * am[i] ≥ 700.0
                    kernel[j,i] = 0.0
                else
                    # de1 is always smaller than 1.0 because τω < 0.0.
                    # de2 could be a large number.
                    de1 = exp(fg[j] * am[i])
                    de2 = exp((fg[j] - β) * am[i])
                    kernel[j,i] = 1.0 / (de1 + de2)
                end
            end
        end
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, fg::FermionicFragmentTimeGrid)

Try to build fermionic kernel function in imaginary time axis. Note that
`fg` contains incomplete imaginary time data.

### Arguments
* am -> Real frequency mesh.
* fg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`FermionicFragmentTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicFragmentTimeGrid)
    ntime = fg.ntime
    nmesh = am.nmesh
    β = fg.β

    kernel = zeros(F64, ntime, nmesh)

    #
    # Old implementation
    # exp(βω ≈ 700.0) will throw an Inf and K becomes NaN.
    # We should avoid this situation.
    #
    # for i = 1:nmesh
    #     de = 1.0 + exp(-β * am[i])
    #     for j = 1:ntime
    #         kernel[j,i] = exp(-fg[j] * am[i]) / de
    #     end
    # end
    #
    # New implementation to avoid numerical overflow
    #
    for i = 1:nmesh
        # ω ≥ 0.0
        if am[i] ≥ 0.0
            de = 1.0 + exp(-β * am[i])
            for j = 1:ntime
                kernel[j,i] = exp(-fg[j] * am[i]) / de
            end
        # ω < 0.0
        else
            for j = 1:ntime
                # Now denominator goes to infinity and K goes to zero.
                if (fg[j] - β) * am[i] ≥ 700.0
                    kernel[j,i] = 0.0
                else
                    # de1 is always smaller than 1.0 because τω < 0.0.
                    # de2 could be a large number.
                    de1 = exp(fg[j] * am[i])
                    de2 = exp((fg[j] - β) * am[i])
                    kernel[j,i] = 1.0 / (de1 + de2)
                end
            end
        end
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)

Try to build fermionic kernel function in Matsubara frequency axis. This
function support the so-called preblur algorithm.

### Arguments
* am -> Real frequency mesh.
* fg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`FermionicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = fg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    # No preblur
    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                _kernel[j,i] = 1.0 / (im * fg[j] - am[i])
            end
        end
    # The preblur trick is used
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
    build_kernel(am::AbstractMesh, fg::FermionicFragmentMatsubaraGrid)

Try to build fermionic kernel function in Matsubara frequency axis. This
function support the so-called preblur algorithm.

### Arguments
* am -> Real frequency mesh.
* fg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, fg::FermionicFragmentMatsubaraGrid)
    blur = get_m("blur")
    nfreq = fg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)

    # No preblur
    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                _kernel[j,i] = 1.0 / (im * fg[j] - am[i])
            end
        end
    # The preblur trick is used
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

### Arguments
* am -> Real frequency mesh.
* bg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        de = 1.0 - exp(-β * am[i])
        for j = 1:ntime
            kernel[j,i] = am[i] * exp(-bg[j] * am[i]) / de
        end
    end
    #
    # Be careful, am[1] is not 0.0!
    # We have to find out where the zero point is in the real mesh.
    _, zero_point = findmin(abs.(am.mesh))
    if am[zero_point] == 0.0
        @. kernel[:,zero_point] = 1.0 / β
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicFragmentTimeGrid)

Try to build bosonic kernel function in imaginary time axis. Note that
`bg` contains incomplete imaginary time data.

### Arguments
* am -> Real frequency mesh.
* bg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicFragmentTimeGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicFragmentTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        de = 1.0 - exp(-β * am[i])
        for j = 1:ntime
            kernel[j,i] = am[i] * exp(-bg[j] * am[i]) / de
        end
    end
    #
    # Be careful, am[1] is not 0.0!
    # We have to find out where the zero point is in the real mesh.
    _, zero_point = findmin(abs.(am.mesh))
    if am[zero_point] == 0.0
        @. kernel[:,zero_point] = 1.0 / β
    end

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis.

### Arguments
* am -> Real frequency mesh.
* bg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    nfreq = bg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)
    #
    for i = 1:nmesh
        for j = 1:nfreq
            _kernel[j,i] = am[i] / (im * bg[j] - am[i])
        end
    end
    #
    # Be careful, am[1] is not 0.0!
    # We have to find out where the zero point is in the real mesh.
    _, zero_point = findmin(abs.(am.mesh))
    if am[zero_point] == 0.0 && bg[1] == 0.0
        _kernel[1,zero_point] = -1.0
    end
    #
    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel(am::AbstractMesh, bg::BosonicFragmentMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis.

### Arguments
* am -> Real frequency mesh.
* bg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function build_kernel(am::AbstractMesh, bg::BosonicFragmentMatsubaraGrid)
    nfreq = bg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)
    #
    for i = 1:nmesh
        for j = 1:nfreq
            _kernel[j,i] = am[i] / (im * bg[j] - am[i])
        end
    end
    #
    # Be careful, am[1] is not 0.0!
    # We have to find out where the zero point is in the real mesh.
    _, zero_point = findmin(abs.(am.mesh))
    if am[zero_point] == 0.0 && bg[1] == 0.0
        _kernel[1,zero_point] = -1.0
    end
    #
    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)

Try to build bosonic kernel function in imaginary time axis (just for
correlator of Hermitian operator only).

### Arguments
* am -> Real frequency mesh.
* bg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicImaginaryTimeGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        r = am[i] / (1.0 - exp(-β * am[i]))
        for j = 1:ntime
            kernel[j,i] = r * (exp(-am[i] * bg[j]) + exp(-am[i] * (β - bg[j])))
        end
    end
    #
    # Perhaps we should check am[1] here!
    @assert am[1] == 0.0
    @. kernel[:,1] = 2.0 / β

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicFragmentTimeGrid)

Try to build bosonic kernel function in imaginary time axis (just for
correlator of Hermitian operator only). Note that `bg` contains
incomplete imaginary time data.

### Arguments
* am -> Real frequency mesh.
* bg -> Imaginary time grid.

### Returns
* kernel -> Kernel function, K(τ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicFragmentTimeGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicFragmentTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    kernel = zeros(F64, ntime, nmesh)
    #
    for i = 1:nmesh
        r = am[i] / (1.0 - exp(-β * am[i]))
        for j = 1:ntime
            kernel[j,i] = r * (exp(-am[i] * bg[j]) + exp(-am[i] * (β - bg[j])))
        end
    end
    #
    # Perhaps we should check am[1] here!
    @assert am[1] == 0.0
    @. kernel[:,1] = 2.0 / β

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis (just
for correlator of Hermitian operator only). This function support the
so-called preblur algorithm.

### Arguments
* am -> Real frequency mesh.
* bg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicMatsubaraGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    blur = get_m("blur")
    nfreq = bg.nfreq
    nmesh = am.nmesh

    kernel = zeros(F64, nfreq, nmesh)

    # No preblur
    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                kernel[j,i] = am[i] ^ 2.0 / ( bg[j] ^ 2.0 + am[i] ^ 2.0 )
                kernel[j,i] = -2.0 * kernel[j,i]
            end
        end
        #
        # Perhaps we should check am[i] and bg[j] here!
        if am[1] == 0.0 && bg[1] == 0.0
            kernel[1,1] = -2.0
        end
    # The preblur trick is used
    else
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)

        I₁ = zeros(F64, nsize)
        I₂ = zeros(F64, nsize)
        I₃ = zeros(F64, nsize)

        for i = 1:nmesh
            for j = 1:nfreq
                g² = bg[j] ^ 2.0
                for k = 1:nsize
                    A² = (bmesh[k] + am[i]) ^ 2.0; I₁[k] = -2.0 * A² / (A² + g²)
                    B² = (bmesh[k] - am[i]) ^ 2.0; I₂[k] = -2.0 * B² / (B² + g²)
                end
                #
                # Perhaps we should check am[i] and bg[j] here!
                if i == 1 && j == 1
                    @assert am[i] == 0.0
                    @assert bg[j] == 0.0
                    I₁ .= -2.0
                    I₂ .= -2.0
                end
                @. I₃ = (I₁ + I₂) * gaussian / 2.0
                kernel[j,i] = simpson(bmesh, I₃)
            end
        end
    end

    return kernel
end

"""
    build_kernel_symm(am::AbstractMesh, bg::BosonicFragmentMatsubaraGrid)

Try to build bosonic kernel function in Matsubara frequency axis (just
for correlator of Hermitian operator only). This function support the
so-called preblur algorithm.

### Arguments
* am -> Real frequency mesh.
* bg -> Matsubara frequency grid.

### Returns
* kernel -> Kernel function, K(iωₙ,ω).

See also: [`AbstractMesh`](@ref), [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function build_kernel_symm(am::AbstractMesh, bg::BosonicFragmentMatsubaraGrid)
    blur = get_m("blur")
    nfreq = bg.nfreq
    nmesh = am.nmesh

    kernel = zeros(F64, nfreq, nmesh)

    # No preblur
    if blur isa Missing || blur < 0.0
        for i = 1:nmesh
            for j = 1:nfreq
                kernel[j,i] = am[i] ^ 2.0 / ( bg[j] ^ 2.0 + am[i] ^ 2.0 )
                kernel[j,i] = -2.0 * kernel[j,i]
            end
        end
        #
        # Perhaps we should check am[i] and bg[j] here!
        if am[1] == 0.0 && bg[1] == 0.0
            kernel[1,1] = -2.0
        end
    # The preblur trick is used
    else
        bmesh, gaussian = make_gauss_peaks(blur)
        nsize = length(bmesh)

        I₁ = zeros(F64, nsize)
        I₂ = zeros(F64, nsize)
        I₃ = zeros(F64, nsize)

        for i = 1:nmesh
            for j = 1:nfreq
                g² = bg[j] ^ 2.0
                for k = 1:nsize
                    A² = (bmesh[k] + am[i]) ^ 2.0; I₁[k] = -2.0 * A² / (A² + g²)
                    B² = (bmesh[k] - am[i]) ^ 2.0; I₂[k] = -2.0 * B² / (B² + g²)
                end
                #
                # Perhaps we should check am[i] and bg[j] here!
                if i == 1 && j == 1
                    @assert am[i] == 0.0
                    @assert bg[j] == 0.0
                    I₁ .= -2.0
                    I₂ .= -2.0
                end
                @. I₃ = (I₁ + I₂) * gaussian / 2.0
                kernel[j,i] = simpson(bmesh, I₃)
            end
        end
    end

    return kernel
end

"""
    make_blur(am::AbstractMesh, A::Vector{F64}, blur::F64)

Try to blur the given spectrum `A`, which is defined in `am`. And `blur`
is the blur parameter.

### Arguments
* am   -> Real frequency mesh.
* A    -> Spectral function.
* blur -> Blur parameter. It must be larger than 0.0.

### Returns
* A    -> It is updated in this function.
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

Perform singular value decomposition for the input matrix `kernel`.

kernel = U Σ Vᵀ

Supposed that kernel is a m × n matrix, then U is m × m, Σ is m × n,
and V is n × n. For Σ, only the diagonal elements are non-zero.

### Arguments
* kernel -> Fermionic or bosonic kernel matrix.

### Returns
* U -> A m × n matrix.
* S -> Diagonal elements of Σ.
* V -> A n × n matrix.
"""
function make_singular_space(kernel::Matrix{F64})
    U, S, V = svd(kernel)

    n_svd = count(x -> x ≥ 1e-10, S)
    U_svd = U[:,1:n_svd]
    V_svd = V[:,1:n_svd]
    S_svd = S[1:n_svd]

    return U_svd, V_svd, S_svd
end

#=
*Remarks* : *About Preblur*

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

Try to generate a gaussian peak along a linear mesh, whose energy range
is `[-5 * blur, +5 * blur]`. The number of mesh points is fixed to 201.

### Arguments
* blur -> This parameter is used to control the width of gaussian peak.

### Returns
* bmesh -> A linear mesh in [-5 * blur, 5 * blur].
* gaussian -> A gaussian peak at `bmesh`.
"""
function make_gauss_peaks(blur::F64)
    @assert blur > 0.0
    nsize = 201
    bmesh = collect(LinRange(-5.0 * blur, 5.0 * blur, nsize))
    norm = 1.0 / (blur * sqrt(2.0 * π))
    gaussian = norm * exp.(-0.5 * (bmesh / blur) .^ 2.0)
    return bmesh, gaussian
end
