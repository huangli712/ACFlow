using LinearAlgebra

𝐴 = [3.8163e-11        0e+00   1.6561e-08;
          0e+00   7.9080e-05        0e+00;
     1.6561e-08        0e+00   7.3933e-06]

𝐵 = [3.8163e-11 0.0 1.6561e-8;
     0.0 7.9080e-5 -0.0;
     1.6561e-8 0.0 7.3933e-6]

𝐶 = [3.8163 0.00 1.6561;
     0.00 7.9080 0.00;
     1.6561 0e+00 7.3933]

𝐷 = [3.8163 0e+00 1.6561;
     0.00 7.9080 -0.00;
     1.6561 0e+00 7.3933]

#@show 𝐴 == 𝐵
𝐹 = eigen(𝐶)
#@show 𝐶
#println(typeof(𝐴))
#println(𝐹.values)
println(𝐹.vectors)
#@show (𝐹.vectors) * diagm(𝐹.values) * (𝐹.vectors)'



F = eigen(𝐷)
#println(typeof(𝐵))
#println(F.values)
println(F.vectors)
#@show 𝐷
#@show (F.vectors) * diagm(F.values) * (F.vectors)'

