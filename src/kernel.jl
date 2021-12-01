#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/01
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

    NCg = 3 * rfg.nul - 1
    Nx = length(rfg.grid)
    B = zeros(F64, NCg, NCg)
    Ps = zeros(F64, NCg, Nx)
    Pg = zeros(F64, 4 * rfg.nul, 4 * rfg.nul - 1)

    B[1,1] = 1.0
	B[1,2] = 1.0
	B[2,1] = 3.0
	B[2,2] = 2.0
	B[2,5] = -RDg[1]
	B[3,1] = 6.0
	B[3,2] = 2.0
	B[3,4] = -2.0 * (RDg[1])^2
	
	Ps[1,1] = -1.0
	Ps[1,2] = 1.0
	
	Pg[1,1] = 1.0
	Pg[2,2] = 1.0
	Pg[4,NCg+1] = 1.0

    for j = 1:Ng-1
        B[3*j+1,3*j+0] = 1.0
        B[3*j+1,3*j+1] = 1.0
        B[3*j+1,3*j+2] = 1.0
        B[3*j+2,3*j+0] = 3.0
        B[3*j+2,3*j+1] = 2.0
        B[3*j+2,3*j+2] = 1.0
        B[3*j+2,3*j+5] = -RDg[j+1]
        B[3*j+3,3*j+0] = 6.0
        B[3*j+3,3*j+1] = 2.0
        B[3*j+3,3*j+4] = -2.0 * (RDg[j+1])^2
            
        Ps[3*j+1,j+1] = -1.0
        Ps[3*j+1,j+2] = 1.0
            
        Pg[4*j+1,3*j+0] = 1.0
        Pg[4*j+2,3*j+1] = 1.0
        Pg[4*j+3,3*j+2] = 1.0
        Pg[4*j+4,NCg+j+1] = 1.0
    end
    
end