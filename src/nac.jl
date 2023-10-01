#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/10/01
#

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

mutable struct NevanlinnaSolver
    Gᵥ::Vector{APC}
    grid::AbstractGrid
    mesh::AbstractMesh
    Gout::Vector{APC}
    phis::Vector{APC}          # phis in schur algorithm
    abcd::Array{APC,3}         # continued fractions
    H_max::I64                 # upper cut off of H
    H_min::I64                 # lower cut off of H
    H::I64                     # current value of H
    ab_coeff::Vector{C64}      # current solution for H
    hardy_matrix::Array{APC,2} # hardy_matrix for H
    iter_tol::I64              # upper bound of iteration
    ini_iter_tol::I64          # upper bound of iteration for H_min
end

function NevanlinnaSolver(
                  wn          ::Vector{APC},
                  gw          ::Vector{APC},
                  N_real      ::I64,
                  H_max       ::I64,
                  iter_tol    ::I64,
                  ;
                  pick_check  ::Bool=true,
                  optimization::Bool=true,
                  ini_iter_tol::I64=500,
                  ham_option  ::Bool=false #option for using in Hamburger moment problem
)
    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(wn) == length(gw)
    N_imag = length(wn) 
    @show N_imag

    if pick_check
        @time opt_N_imag =  calc_Nopt(wn, gw)
    else 
        opt_N_imag = N_imag
    end
    @show opt_N_imag

    β::APF = 100.0
    grid = FermionicMatsubaraGrid(opt_N_imag, β, reverse(imag.(wn[1:opt_N_imag])))
    mesh = make_mesh(T = APF)
    Gout = zeros(APC, N_real)
    Gᵥ = calc_mobius(-gw[1:opt_N_imag])
    reverse!(Gᵥ)

    phis = calc_phis(grid, Gᵥ)
    abcd = calc_abcd(grid, mesh, phis)

    H_min::Int64 = 1
    ab_coeff = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(mesh, H_min)

    sol = NevanlinnaSolver(Gᵥ, grid, mesh, Gout, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, ini_iter_tol)

    if ham_option
        return sol
    end
    
    if optimization
        calc_H_min(sol)
    else
        evaluation!(sol)
    end

    return sol
end

function calc_mobius(z::Vector{APC})
    _z = similar(z)
    @. _z = (z - im) / (z + im)
    return _z
end

function calc_pick(k::I64, λ::Vector{APC}, ℎ::Vector{APC})
    pick = zeros(APC, k, k)
    for j = 1:k
        for i = 1:k
            num = one(APC) - λ[i] * conj(λ[j])
            den = one(APC) - ℎ[i] * conj(ℎ[j])
            pick[i,j] = num / den
        end
        pick[j,j] += APC(1e-250)
    end
    return issuccess(cholesky(pick, check = false))
end

function calc_Nopt(wn::Vector{APC}, gw::Vector{APC})
    N = length(wn)

    freq = calc_mobius(wn)
    val = calc_mobius(-gw)

    k = 0
    success = true
    while success && k ≤ N
        k += 1
        success = calc_pick(k, val, freq)
    end

    if !success
        println("N_imag is setted as $(k-1)")
        return k-1
    else
        println("N_imag is setted as $(N)")
        return N
    end
end

function calc_phis(grid::AbstractGrid, Gᵥ::Vector{APC})
    Nopt = length(grid)

    Φ = zeros(APC, Nopt) 
    𝒜 = zeros(APC, 2, 2, Nopt)
    ∏ = zeros(APC, 2, 2)
    𝑔 = grid.ω * im

    for i in 1:Nopt
        view(𝒜,:,:,i) .= Matrix{APC}(I, 2, 2)
    end

    Φ[1] = Gᵥ[1]
    for j in 1:Nopt-1
        for k in j+1:Nopt
            ∏[1,1] = ( 𝑔[k] - 𝑔[j] ) / ( 𝑔[k] - conj(𝑔[j]) )
            ∏[1,2] = Φ[j]
            ∏[2,1] = conj(Φ[j]) * ( 𝑔[k] - 𝑔[j] ) / ( 𝑔[k] - conj(𝑔[j]) )
            ∏[2,2] = one(APC)
            view(𝒜,:,:,k) .= view(𝒜,:,:,k) * ∏
        end
        num = 𝒜[1,2,j+1] - 𝒜[2,2,j+1] * Gᵥ[j+1]
        den = 𝒜[2,1,j+1] * Gᵥ[j+1] - 𝒜[1,1,j+1]
        Φ[j+1] = num / den
    end

    return Φ
end

function calc_abcd(grid::AbstractGrid, mesh::AbstractMesh, Φ::Vector{APC})
    eta::APF = get_n("eta")

    ngrid = length(grid)
    nmesh = length(mesh)

    𝒜 = zeros(APC, 2, 2, nmesh)
    𝑔 = grid.ω * im
    𝑚 = mesh.mesh .+ im * eta

    for i in 1:nmesh
        result = Matrix{APC}(I, 2, 2)
        𝑧 = 𝑚[i]
        for j in 1:ngrid
            ∏ = zeros(APC, 2, 2)
            ∏[1,1] = ( 𝑧 - 𝑔[j] ) / ( 𝑧 - conj(𝑔[j]) )
            ∏[1,2] = Φ[j]
            ∏[2,1] = conj(Φ[j]) * ( 𝑧 - 𝑔[j] ) / ( 𝑧 - conj(𝑔[j]) )
            ∏[2,2] = one(APC)
            result *= ∏
        end

        𝒜[:,:,i] .= result
    end

    return 𝒜
end

function check_causality(hardy_matrix::Array{APC,2}, ab_coeff::Vector{C64})
    param = hardy_matrix*ab_coeff

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

function evaluation!(sol::NevanlinnaSolver)
    causality = check_causality(sol.hardy_matrix, sol.ab_coeff)
    if causality
        param = sol.hardy_matrix*sol.ab_coeff
        theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
        sol.Gout .= im * (one(APC) .+ theta) ./ (one(APC) .- theta)
    end

    return causality
end

function hardy_basis(z::APC, k::I64)
    w = (z-im)/(z+im)
    0.5*im*(w^(k+1)-w^k)/(sqrt(pi))
end

function calc_hardy_matrix(am::AbstractMesh, H::I64)
    N_real = length(am)
    hardy_matrix = zeros(APC, N_real, 2*H)
    eta::APF = get_n("eta")
    freq = am.mesh .+ eta * im
    for k in 1:H
        hardy_matrix[:,2*k-1] .=      hardy_basis.(freq,k-1)
        hardy_matrix[:,2*k]   .= conj(hardy_basis.(freq,k-1))
    end
    return hardy_matrix
end

function calc_H_min(sol::NevanlinnaSolver,)::Nothing
    H_bound::Int64 = 50
    for iH in 1:H_bound
        println("H=$(iH)")
        zero_ab_coeff = zeros(C64, 2*iH)

        causality, optim = hardy_optim!(sol, iH, zero_ab_coeff, iter_tol=sol.ini_iter_tol)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            sol.H_min = sol.H
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end

        if iH == H_bound
            error("H_min does not exist")
        end
    end
end

function calc_functional(sol::NevanlinnaSolver, H::Int64, ab_coeff::Vector{C64}, hardy_matrix::Array{APC,2})
    param = hardy_matrix*ab_coeff

    theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
    green = im * (one(APC) .+ theta) ./ (one(APC) .- theta)
    A = F64.(imag(green)./pi)

    tot_int = trapz(sol.mesh, A)
    second_der = integrate_squared_second_deriv(sol.mesh.mesh, A) 

    alpha = get_n("alpha")
    func = abs(1.0-tot_int)^2 + alpha*second_der

    return func
end

function hardy_optim!(
                sol::NevanlinnaSolver,
                H::I64,
                ab_coeff::Vector{C64};
                iter_tol::I64=sol.iter_tol,
                )::Tuple{Bool, Bool}
    loc_hardy_matrix = calc_hardy_matrix(sol.mesh, H)

    function functional(x::Vector{C64})::F64
        return calc_functional(sol, H, x, loc_hardy_matrix)
    end

    function jacobian(J::Vector{C64}, x::Vector{C64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = true))
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(loc_hardy_matrix, Optim.minimizer(res))

    if causality && (Optim.converged(res))
        sol.H = H
        sol.ab_coeff = Optim.minimizer(res)
        sol.hardy_matrix = loc_hardy_matrix
        evaluation!(sol)
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

function solve(S::NevanACSolver, rd::RawData)
    N_real    = 1000  #demension of array of output
    H_max     = 12    #cutoff of Hardy basis
    iter_tol  = 1000  #upper bound of iteration

    T = BigFloat
    setprecision(128)

    input_smpl = zeros(Complex{BigFloat},52)
    input_gw = zeros(Complex{BigFloat},52)

    dlm = readdlm("gw.data")
    @show size(dlm), typeof(dlm)
    @. input_smpl = dlm[:,1] * im
    @. input_gw = dlm[:,2] + dlm[:,3] * im
    wo_sol = NevanlinnaSolver(input_smpl, input_gw, N_real, H_max, iter_tol)

    open("twopeak_wo_opt.dat","w") do f
        for i in 1:N_real
            println(f, "$(F64(wo_sol.mesh[i]))",  "\t", "$(F64(imag.(wo_sol.Gout[i]/pi)))")
        end
    end
end
