#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/30
#

function calc_kernel(ω::FermionicMatsubaraGrid, rfg::RealFrequencyGrid)
    println("here in calc_kernel()")
    @show length(ω.grid)
    @show length(rfg.grid)

    dul = rfg.dw / (rfg.wl - rfg.dw - rfg.w0l) / (rfg.wl - rfg.w0l)
    ul = F64[]
    for i = 1:rfg.nul
        push!(ul, -dul * i)
    end

    dur = rfg.dw / (rfg.wr + rfg.dw - rfg.w0r) / (rfg.wr - rfg.w0r)
    ur = F64[]
    for i = 1:rfg.nur
        insert!(ur, 1, dur * i)
    end

    #@show ur
    #@show rfg.w0l, rfg.w0r

    spline_matrix(rfg)
end

function spline_matrix(rfg::RealFrequencyGrid)
    println("here is spline_matrix")

    v1 = rfg.grid[3:rfg.nul+1] .- rfg.w0l
    v2 = rfg.grid[1:rfg.nul-1] .- rfg.w0l
    v3 = rfg.grid[2:rfg.nul] .- rfg.grid[1:rfg.nul-1]
    v4 = rfg.grid[3:rfg.nul+1] .- rfg.grid[2:rfg.nul]

    RDg = (v1 ./ v2) .* (v3 ./ v4)
    #@show RDg

    
end