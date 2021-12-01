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
    Mg = spline_matrix_g(rfg)
    #@show Mg
    spline_matrix_c(rfg)
end

function spline_matrix_g(rfg::RealFrequencyGrid)
    # println("here is spline_matrix")

    v1 = rfg.grid[3:rfg.nul+1] .- rfg.w0l
    v2 = rfg.grid[1:rfg.nul-1] .- rfg.w0l
    v3 = rfg.grid[2:rfg.nul] .- rfg.grid[1:rfg.nul-1]
    v4 = rfg.grid[3:rfg.nul+1] .- rfg.grid[2:rfg.nul]

    RDg = (v1 ./ v2) .* (v3 ./ v4)
    #@show RDg

    Ng = rfg.nul
    NCg = 3 * Ng - 1
    Nx = length(rfg.grid)
    B = zeros(F64, NCg, NCg)
    Ps = zeros(F64, NCg, Nx)
    Pg = zeros(F64, 4 * Ng, 4 * Ng - 1)

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

    for j = 1:Ng-2
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

    j = Ng - 1
    B[3*j+1,3*j+0] = 1.0
	B[3*j+1,3*j+1] = 1.0
	B[3*j+1,3*j+2] = 1.0
	B[3*j+2,3*j+0] = 3.0
	B[3*j+2,3*j+1] = 2.0
	B[3*j+2,3*j+2] = 1.0

    fdAg = ( rfg.grid[j+2] - rfg.w0l ) / ( rfg.grid[j+1] - rfg.w0l )
    fdAg = ( rfg.grid[j+2] - rfg.grid[j+1] ) / ( rfg.grid[j+3] - rfg.grid[j+1] ) * fdAg
    @show fdAg

	Ps[3*j+1,j+1] = -1.0
	Ps[3*j+1,j+2] = 1.0
	Ps[3*j+2,j+1] = -fdAg
	Ps[3*j+2,j+3] = fdAg
	
	Pg[4*j+1,3*j+0] = 1.0
	Pg[4*j+2,3*j+1] = 1.0
	Pg[4*j+3,3*j+2] = 1.0
	Pg[4*j+4,NCg+j+1] = 1.0

    IB = Matrix{F64}(I, NCg, NCg)
    invB = B \ IB
    #@show invB

    IA = Matrix{F64}(I, Nx, Nx)
    PA = IA[1:Ng,1:Nx]
    Lg = vcat(invB * Ps, PA)
    #@show Lg
    Mg = Pg * Lg
    #@show size(Mg)
    #@show Mg
    return Mg
end

function spline_matrix_c(rfg::RealFrequencyGrid)
    nuc = length(rfg.grid) - rfg.nur - rfg.nul
    @show rfg.nul, nuc, nuc + rfg.nul
    v1 = rfg.grid[rfg.nul + 2 : rfg.nul + nuc - 1]
    v2 = rfg.grid[rfg.nul + 1 : rfg.nul + nuc - 2]
    v3 = rfg.grid[rfg.nul + 3 : rfg.nul + nuc - 0]
    v4 = rfg.grid[rfg.nul + 2 : rfg.nul + nuc - 1]
    #@show v4

    RDc = (v1 .- v2) ./ (v3 .- v4)
    #@show RDc

    Ng = rfg.nul
    Nc = nuc - 1
    NCc = 3 * Nc - 1
    Nx = length(rfg.grid)
    @show NCc, Nc, Nx

    B = zeros(F64, NCc, NCc)
    Ps = zeros(F64, NCc, Nx)
    Pg = zeros(F64, 4 * Nc, 4 * Nc)

    B[1,1] = 1.0
	B[1,2] = 1.0
	B[2,1] = 3.0
	B[2,2] = 2.0
	B[2,5] = -RDc(0)
	B[3,1] = 6.0
	B[3,2] = 2.0
	B[3,4] = -2.0 * pow(RDc(0),2)
	
	Ps[1,Ng+0] = ( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[1,Ng+1] = -1.0
	Ps[1,Ng+2] = 1.0 - ( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[2,Ng+0] = +( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[2,Ng+2] = -( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	
	Pg[1,1] = 1.0
	Pg[2,2] = 1.0
	Pg[3,NCc+1] = -( rfg.grid[Ng+2] - rfg.grid[Ng+1]) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Pg[3,NCc+3] = +( rfg.grid[Ng+2] - rfg.grid[Ng+1]) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Pg[4,NCc+2] = 1.0

end

function spline_matrix_d(rfg::RealFrequencyGrid)
end