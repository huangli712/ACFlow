
function make_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)
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

#=
```math
\int
```
=#

function make_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)
    nfreq = fg.nfreq
    nmesh = am.nmesh

    _kernel = zeros(C64, nfreq, nmesh)
    for i = 1:nmesh
        for j = 1:nfreq
            _kernel[j,i] = 1.0 / (im * fg[j] - am[i])
        end
    end

    kernel = vcat(real(_kernel), imag(_kernel))

    return kernel
end

function make_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
    ntime = bg.ntime
    nmesh = am.nmesh
    β = bg.β

    #@show am.mesh
    #@show bg.τ
    #error()
    kernel = zeros(F64, ntime, nmesh)
    for i = 1:nmesh
        r = 0.25 * am[i] / (1. - exp(-β * am[i]))
        #@show r
        for j = 1:ntime
            kernel[j,i] = r * (exp(-am[i] * bg[j]) + exp(-am[i] * (β - bg[j]))) 
        end
    end

    @. kernel[:,1] = 1.0
    return kernel
end

function make_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)
    nfreq = bg.nfreq
    nmesh = am.nmesh

    kernel = zeros(F64, nfreq, nmesh)
    for i = 1:nmesh
        for j = 1:nfreq
            kernel[j,i] = am[i] ^ 2.0 / ( bg[j] ^ 2.0 + am[i] ^ 2.0 )
        end
    end

    if am[1] == 0.0 && bg[1] == 0.0
        kernel[1,1] = 1.0
    end

    return kernel
end

function make_singular_space(kernel::Matrix{F64})
    U, S, V = svd(kernel)
    n_svd = count(x -> x ≥ 1e-10, S)
    U_svd = U[:,1:n_svd]
    V_svd = V[:,1:n_svd]
    S_svd = S[1:n_svd]

    return U_svd, V_svd, S_svd
end