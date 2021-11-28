#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/28
#

function calc_mesh(𝑀::MomentsData)
    SC = 𝑀.𝑀₁ / 𝑀.𝑀₀
    SW = sqrt(𝑀.𝑀₂ / 𝑀.𝑀₀ - SC^2) * 3
    wl = SC - SW / 2
    wr = SC + SW / 2
    SW = wr - wl
    println("SC: ", SC)
    println("SW: ", SW)
    println("wl: ", wl)
    println("wr: ", wr)
end
