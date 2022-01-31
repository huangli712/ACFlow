#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/01/29
#

function blur_kernel(kernel::Matrix{F64}, grid::AbstractGrid, mesh::AbstractMesh)
    blur_width = 0.3
    nsize = 201
    w_int = collect(LinRange(-5.0 * blur_width, 5.0 * blur_width, nsize))
    #@show w_int
    norm = 1.0 / (blur_width * sqrt(2.0 * π))

    #@show norm
    gaussian = norm * exp.(-0.5 * (w_int / blur_width) .^ 2.0)
    #@show gaussian
    #error()

    Mg = reshape(gaussian, (1, 1, nsize))
    MG = reshape(grid.ω, (grid.nfreq, 1, 1))
    Mm = reshape(mesh.mesh, (1, mesh.nmesh, 1))
    Mw = reshape(w_int, (1, 1, nsize))

    integrand = Mg ./ (im * MG .- Mm .- Mw)

    #@show size(integrand), size(kernel)
    #@show grid.nfreq, mesh.nmesh
    new_kernel = zeros(C64, grid.nfreq, mesh.nmesh)
    for j = 1:mesh.nmesh
        for i = 1:grid.nfreq
            new_kernel[i,j] = simpson(w_int, integrand[i,j,:])
            @show i, j, new_kernel[i,j]
        end
    end
    error()
end

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