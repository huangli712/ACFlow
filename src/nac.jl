#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/09/27
#

push!(LOAD_PATH, "/home/lihuang/Downloads/Nevanlinna.jl-main/src")
using Nevanlinna

function solve(S::NevanACSolver, rd::RawData)
    N_real    = 1000  #demension of array of output
    omega_max = 10.0  #energy cutoff of real axis
    eta       = 0.001 #broaden parameter 
    sum_rule  = 1.0   #sum rule
    H_max     = 12    #cutoff of Hardy basis
    lambda    = 1e-4  #regularization parameter
    iter_tol  = 1000  #upper bound of iteration

    T = BigFloat
    setprecision(128)

    input_smpl = zeros(Complex{BigFloat},52)
    input_gw = zeros(Complex{BigFloat},52)

    dlm = readdlm("gw.data")
    @show size(dlm)
    exit(-1)
    wo_sol = Nevanlinna.NevanlinnaSolver(input_smpl, input_gw, N_real, omega_max, eta, sum_rule, H_max, iter_tol, lambda, verbose=true)
    exit(-1)
end