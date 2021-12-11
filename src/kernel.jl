#
# Project : Gardenia
# Source  : kernel.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/10
#

function calc_kernel(ω::FermionicMatsubaraGrid, rfg::RealFrequencyGrid)
    #println("here in calc_kernel()")
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

    MM = _spline_matrix(rfg)
    #println("here")

    _kernel_p_g(rfg)
    _kernel_p_c(rfg)
    _kernel_p_d(rfg)

    _kernel_k_g(ul, ω, rfg)
    _kernel_k_c(ω, rfg)
    _kernel_k_d(ur, ω, rfg)
end

function calc_kernel_p()
end

function calc_kernel_k()
end

function calc_kernel_m()
end

function _kernel_p_g(rfg::RealFrequencyGrid)
    #println("here")
    Nw = length(rfg.grid)
    NCfs = 4 * (Nw - 1)
    #@show NCfs
    Nintg = rfg.nul
    #@show Nintg

    Pa_g = zeros(F64, Nintg, NCfs)
    Pb_g_r = zeros(F64, Nintg, NCfs)
    Pc_g_r = zeros(F64, Nintg, NCfs)
    Pd_g_r = zeros(F64, Nintg, NCfs)

    for j = 0:Nintg-1
        Pa_g[j+1,4*j+1] = 1.0
        Pb_g_r[j+1,4*j+2] = 1.0
        Pc_g_r[j+1,4*j+3] = 1.0
        Pd_g_r[j+1,4*j+4] = 1.0
    end
    
    vDG = ( rfg.grid[2:rfg.nul+1] .- rfg.w0l ) .* ( rfg.grid[1:rfg.nul+0] .- rfg.w0l )
    vDG = vDG ./ ( rfg.grid[1:rfg.nul+0] .- rfg.grid[2:rfg.nul+1] )
    #@show vDG
    #@show length(vDG)
    DG = diagm(vDG)

    vDU = rfg.grid[2:rfg.nul+1] .- rfg.w0l
    vDU = vDU ./ ( rfg.grid[1:rfg.nul+0] - rfg.grid[2:rfg.nul+1] )
    #@show vDU
    #@show length(vDU)
    DU = diagm(vDU)

    Pb_g = Pb_g_r - 3.0 * DU * Pa_g
    #@show Pb_g

    Pc_g = Pc_g_r + 3.0 * (DU ^ 2.0) * Pa_g - 2.0 * DU * Pb_g_r
    #@show Pc_g

    Pd_g = Pd_g_r - (DU ^ 3.0) * Pa_g + (DU ^ 2.0) * Pb_g_r - DU * Pc_g_r
    #@show Pd_g

    Pa_g = (DG ^ 3.0) * Pa_g
    Pb_g = (DG ^ 2.0) * Pb_g
    Pc_g = DG * Pc_g
    #@show Pc_g
end

function _kernel_p_c(rfg::RealFrequencyGrid)
    #println("in _kernel_p_c")
    Nw = length(rfg.grid)
    NCfs = 4 * (Nw - 1)
    Nwc = Nw - rfg.nur - rfg.nul
    Nintc = Nwc - 1
    NCg = 4 * rfg.nul
    #@show Nintc, NCfs

    Pa_c = zeros(F64, Nintc, NCfs)
    Pb_c = zeros(F64, Nintc, NCfs)
    Pc_c = zeros(F64, Nintc, NCfs)
    Pd_c = zeros(F64, Nintc, NCfs)

    for j = 0:Nintc-1
        Pa_c[j+1,4*j+1+NCg] = 1.0
		Pb_c[j+1,4*j+2+NCg] = 1.0
		Pc_c[j+1,4*j+3+NCg] = 1.0
		Pd_c[j+1,4*j+4+NCg] = 1.0
    end

    vDC = 1.0 ./ (rfg.grid[rfg.nul+2:Nwc+rfg.nul+0] .- rfg.grid[rfg.nul+1:Nwc+rfg.nul-1])
    #@show vDC
    #@show length(vDC)

    DC = diagm(vDC)

    Pa_c = (DC ^ 3.0) * Pa_c
    Pb_c = (DC ^ 2.0) * Pb_c
    Pc_c = DC * Pc_c

    #@show Pd_c
end

function _kernel_p_d(rfg::RealFrequencyGrid)
    #println("in _kernel_p_d()")
    Nw = length(rfg.grid)
    NCfs = 4 * (Nw - 1)
    Nintd = rfg.nur
    #@show Nintd, NCfs
    Nwc = Nw - rfg.nur - rfg.nul
    Nintc = Nwc - 1
    NCg = 4 * rfg.nul
    NCgc = NCg + 4 * Nintc

    #@show NCgc

    Pa_d = zeros(F64, Nintd, NCfs)
    Pb_d_r = zeros(F64, Nintd, NCfs)
    Pc_d_r = zeros(F64, Nintd, NCfs)
    Pd_d_r = zeros(F64, Nintd, NCfs)

    for j = 0:Nintd-1
		Pa_d[j+1,4*j+1+NCgc] = 1.0
		Pb_d_r[j+1,4*j+2+NCgc] = 1.0
		Pc_d_r[j+1,4*j+3+NCgc] = 1.0
		Pd_d_r[j+1,4*j+4+NCgc] = 1.0
    end

    #Nwc + rfg.nul - 1
    vDD = (rfg.grid[Nwc+rfg.nul+1:Nw+0] .- rfg.w0r) .* (rfg.grid[Nwc+rfg.nul+0:Nw-1] .- rfg.w0r)
    vDD = vDD ./ ( rfg.grid[Nwc+rfg.nul+1:Nw+0] - rfg.grid[Nwc+rfg.nul+0:Nw-1] )
    #@show vDD
    DD = diagm(vDD)

    vDU = rfg.grid[Nwc+rfg.nul+0:Nw-1] .- rfg.w0r
    vDU = vDU ./ ( rfg.grid[Nwc+rfg.nul+1:Nw-0] .- rfg.grid[Nwc+rfg.nul+0:Nw-1] )
    #@show vDU
    DU = diagm(vDU)

    Pb_d = Pb_d_r - 3.0 * DU * Pa_d
    Pc_d = Pc_d_r + 3.0 * (DU ^ 2.0) * Pa_d - 2.0 * DU * Pb_d_r
    Pd_d = Pd_d_r - (DU ^ 3.0) * Pa_d + (DU ^ 2.0) * Pb_d_r - DU * Pc_d_r
    
    Pa_d = (DD ^ 3.0) * Pa_d
    Pb_d = (DD ^ 2.0) * Pb_d
    Pc_d = DD * Pc_d
    #@show Pd_d
end

function _kernel_k_g(ug, ω::FermionicMatsubaraGrid, rfg::RealFrequencyGrid)
    Nn = length(ω.grid)
    Nintg = rfg.nul
    @show Nn, Nintg

    #Ka_g = zeros(C64, Nn, Nintg)
    #Kb_g = zeros(C64, Nn, Nintg)
    #Kc_g = zeros(C64, Nn, Nintg)
    #Kd_g = zeros(C64, Nn, Nintg)

    ug2 = copy(ug)
    push!(ug2, 1.0 / (rfg.wl - rfg.w0l))
    @show length(ug2)

    Wng = zeros(F64, Nn, Nintg)
    Ug = zeros(F64, Nn, length(ug2))

    for i = 1:Nn
        for j = 1:Nintg
            Wng[i,j] = ω.grid[i]
        end
        for j = 1:length(ug2)
            Ug[i,j] = ug2[j]
        end
    end
    #@show size(Ug)
    #@show Ug
    
    #mat atang=atan(
    #    (Wng % ( Ug.cols(1,Nintg) - Ug.cols(0,Nintg-1) ))
        #@show Wng .* (Ug[:,2:Nintg+1] - Ug[:,1:Nintg+0])
    #    /
    #    ( 1 + w0l * ( Ug.cols(1,Nintg) + Ug.cols(0,Nintg-1) ) + ( pow(w0l,2) + pow(Wng,2) ) % Ug.cols(0,Nintg-1) % Ug.cols(1,Nintg) )
    atang = 1.0 .+ rfg.w0l .* ( Ug[:,2:Nintg+1] + Ug[:,1:Nintg+0] )
    atang = atang + ( (rfg.w0l ^ 2.0) .+ (Wng .^ 2.0) ) .* Ug[:,1:Nintg+0] .* Ug[:,2:Nintg+1]
    atang = Wng .* (Ug[:,2:Nintg+1] - Ug[:,1:Nintg+0]) ./ atang
    atang = atan.(atang)
    #@show atang
    #)

    #logg = log(
    logg1 = 1.0 .+ 2.0 * rfg.w0l .* Ug[:,2:Nintg+1] .+ ((Ug[:,2:Nintg+1]) .^ 2.0) .* ((Wng .^ 2.0) .+ (rfg.w0l)^2.0)
    logg2 = 1.0 .+ 2.0 * rfg.w0l .* Ug[:,1:Nintg+0] .+ ((Ug[:,1:Nintg+0]) .^ 2.0) .* ((Wng .^ 2.0) .+ (rfg.w0l)^2.0)
    logg = log.(logg1 ./ logg2)
    #);
    #@show logg

    Ka_g = -( Ug[:,2:Nintg+1] - Ug[:,1:Nintg+0] ) ./ ((Wng .+ im * rfg.w0l) .^ 2.0)
    #@show Ka_g
    Ka_g = Ka_g .- im .* ( ((Ug[:,2:Nintg+1]) .^ 2.0) - ((Ug[:,1:Nintg+0]) .^ 2.0) ) ./ ( 2.0 .* ( Wng .+ im * rfg.w0l ) )
    #@show Ka_g
    Ka_g = Ka_g .+ atang ./ ((Wng .+ im * rfg.w0l) .^ 3.0)
    Ka_g = Ka_g .+ im .* logg ./ (2.0 .* ((Wng .+ im * rfg.w0l) .^ 3.0))
    Ka_g = -Ka_g / (2.0 * π)
    #@show Ka_g

    Kb_g = -im .* ( Ug[:,2:Nintg+1] .- Ug[:,1:Nintg+0] ) ./ (Wng .+ im * rfg.w0l)
    Kb_g = Kb_g .+ im .* atang ./ ((Wng .+ im * rfg.w0l) .^ 2.0)
    Kb_g = Kb_g .- logg ./ ( 2.0 .* ((Wng .+ im * rfg.w0l) .^ 2.0) )
	Kb_g = -Kb_g / (2.0 * π)
    #@show Kb_g

    Kc_g = -atang ./ (Wng .+ im * rfg.w0l) .- im .* logg ./ ( 2.0 .* (Wng .+ im * rfg.w0l) )
    Kc_g = -Kc_g / (2.0 * π)
    #@show Kc_g
	
	Kd_g = -im .* atang .+ logg ./ 2 .- log.(Ug[:,2:Nintg+1] ./ Ug[:,1:Nintg+0])
    Kd_g = -Kd_g / (2.0 * π)
    #@show Kd_g
end

function _kernel_k_c(ω::FermionicMatsubaraGrid, rfg::RealFrequencyGrid)
    Nn = length(ω.grid)
    Nintc = length(rfg.grid) - rfg.nur - rfg.nul - 1
    @show Nn, Nintc

    Ka_c = zeros(C64, Nn, Nintc)
    Kb_c = zeros(C64, Nn, Nintc)
    Kc_c = zeros(C64, Nn, Nintc)
    Kd_c = zeros(C64, Nn, Nintc)

    wc = rfg.grid[rfg.nul+1:rfg.nul+Nintc+1]
    Wnc = zeros(F64, Nn, Nintc)
    Wc = zeros(F64, Nn, length(wc))

    for i = 1:Nn
        for j = 1:Nintc
            Wnc[i,j] = ω.grid[i]
        end
        for j = 1:length(wc)
            Wc[i,j] = wc[j]
        end
    end
    #@show Wnc
    #@show wc
    #@show Wc

    logc = (Wnc .^ 2.0) .+ ((Wc[:,2:Nintc+1]) .^ 2.0)
    logc = logc ./ ((Wnc .^ 2.0) + ((Wc[:,1:Nintc+0]) .^ 2.0))
    logc = log.(logc)
    #@show logc

    atanc2 = Wnc .* (Wc[:,2:Nintc+1] .- Wc[:,1:Nintc+0])
    atanc2 = atanc2 ./ (Wc[:,2:Nintc+1] .* Wc[:,1:Nintc+0] .+ (Wnc .^ 2.0))
    atanc2 = atan.(atanc2)
    #@show atanc2

    logc2 = logc ./ 2.0 .+ im .* atanc2
    dWn = im .* Wnc .- Wc[:,1:Nintc+0]
    dWc = Wc[:,2:Nintc+1] .- Wc[:,1:Nintc+0]

    Ka_c = -(dWn .^ 2.0) .* dWc 
    Ka_c = Ka_c .- dWn .* (dWc .^ 2.0) ./ 2.0
    Ka_c = Ka_c .- (dWc .^ 3.0) ./ 3.0
    Ka_c = Ka_c .- (dWn .^ 3.0) .* logc2

    Kb_c = -dWn .* dWc - (dWc .^ 2.0) ./ 2.0 .- (dWn .^ 2.0) .* logc2
	Kc_c = -dWc .- dWn .* logc2
	Kd_c = -logc2

    Ka_c = Ka_c ./ (2.0 * π)
    Kb_c = Kb_c ./ (2.0 * π)
    Kc_c = Kc_c ./ (2.0 * π)
    Kd_c = Kd_c ./ (2.0 * π)
    #@show Kd_c
end

function _kernel_k_d(ud, ω::FermionicMatsubaraGrid, rfg::RealFrequencyGrid)
end

function _kernel_m_g()
end

function _kernel_m_c()
end

function _kernel_m_d()
end

function _spline_matrix(rfg::RealFrequencyGrid)
    Mg = _spline_matrix_g(rfg)
    #@show Mg
    Mc = _spline_matrix_c(rfg)
    Md = _spline_matrix_d(rfg)

    MM = vcat(Mg, Mc, Md)
    #@show size(MM)
    #@show MM
    return MM
end

function _spline_matrix_g(rfg::RealFrequencyGrid)
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
    #@show fdAg

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

function _spline_matrix_c(rfg::RealFrequencyGrid)
    nuc = length(rfg.grid) - rfg.nur - rfg.nul
    #@show rfg.nul, nuc, nuc + rfg.nul
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
    #@show NCc, Nc, Nx

    B = zeros(F64, NCc, NCc)
    Ps = zeros(F64, NCc, Nx)
    Pc = zeros(F64, 4 * Nc, 4 * Nc)

    B[1,1] = 1.0
	B[1,2] = 1.0
	B[2,1] = 3.0
	B[2,2] = 2.0
	B[2,5] = -RDc[1]
	B[3,1] = 6.0
	B[3,2] = 2.0
	B[3,4] = -2.0 * (RDc[1])^2
	
	Ps[1,Ng+0] = ( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[1,Ng+1] = -1.0
	Ps[1,Ng+2] = 1.0 - ( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[2,Ng+0] = +( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Ps[2,Ng+2] = -( rfg.grid[Ng+2] - rfg.grid[Ng+1] ) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	
	Pc[1,1] = 1.0
	Pc[2,2] = 1.0
	Pc[3,NCc+1] = -( rfg.grid[Ng+2] - rfg.grid[Ng+1]) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Pc[3,NCc+3] = +( rfg.grid[Ng+2] - rfg.grid[Ng+1]) / ( rfg.grid[Ng+2] - rfg.grid[Ng] )
	Pc[4,NCc+2] = 1.0

	for j = 1:Nc-2
        B[3*j+1,3*j+0] = 1.0
        B[3*j+1,3*j+1] = 1.0
        B[3*j+1,3*j+2] = 1.0
        B[3*j+2,3*j+0] = 3.0
        B[3*j+2,3*j+1] = 2.0
        B[3*j+2,3*j+2] = 1.0
        B[3*j+2,3*j+5] = -RDc[j+1]
        B[3*j+3,3*j+0] = 6.0
        B[3*j+3,3*j+1] = 2.0
        B[3*j+3,3*j+4] = -2.0 * (RDc[j+1])^2
            
        Ps[3*j+1,Ng+j+1] = -1.0
        Ps[3*j+1,Ng+j+2] = +1.0
            
        Pc[4*j+1,3*j+0] = 1.0
        Pc[4*j+2,3*j+1] = 1.0
        Pc[4*j+3,3*j+2] = 1.0
        Pc[4*j+4,NCc+j+2] = 1.0
    end
    
	j = Nc - 1
	B[3*j+1,3*j+0] = 1
	B[3*j+1,3*j+1] = 1
	B[3*j+1,3*j+2] = 1
	B[3*j+2,3*j+0] = 3
	B[3*j+2,3*j+1] = 2
	B[3*j+2,3*j+2] = 1
	
    fdAc = rfg.grid[Ng + Nc + 1] - rfg.grid[Ng + Nc]
    fdAc = fdAc / ( rfg.grid[Ng + Nc + 2] - rfg.grid[Ng + Nc] )
    #@show fdAc

	Ps[3*j+1,Ng+j+1] = -1
	Ps[3*j+1,Ng+j+2] = 1
	Ps[3*j+2,Ng+j+1] = -fdAc
	Ps[3*j+2,Ng+j+3] = +fdAc
	
	Pc[4*j+1,3*j+0] = 1.0
	Pc[4*j+2,3*j+1] = 1.0
	Pc[4*j+3,3*j+2] = 1.0
	Pc[4*j+4,NCc+j+2] = 1.0

    IB = Matrix{F64}(I, NCc, NCc)
    invB = B \ IB
    #@show invB

    IA = Matrix{F64}(I, Nx, Nx)
    PA = IA[Ng:Ng+Nc,1:Nx]
    Lc = vcat(invB * Ps, PA)
    #@show Lc
    Mc = Pc * Lc
    #@show size(Mc)
    #@show Mc
    return Mc
end

function _spline_matrix_d(rfg::RealFrequencyGrid)
    Nx = length(rfg.grid)
    Nd = rfg.nur
    NCd = 3 * Nd - 1

    #Nx - Nd - 1
    v1 = rfg.grid[Nx - Nd + 2 : Nx - 0] .- rfg.grid[Nx - Nd + 1 : Nx - 1]
    v2 = rfg.grid[Nx - Nd + 1 : Nx - 1] .- rfg.grid[Nx - Nd + 0 : Nx - 2]
    v3 = rfg.grid[Nx - Nd + 0 : Nx - 2] .- rfg.w0r
    v4 = rfg.grid[Nx - Nd + 2 : Nx - 0] .- rfg.w0r
    RDd = (v1 ./ v2) .* (v3 ./ v4)
    #@show RDd
    #@show length(RDd)

    fdAd = ( rfg.grid[Nx - Nd] - rfg.w0r ) / ( rfg.grid[Nx - Nd + 1] - rfg.w0r )
    fdAd = fdAd * ( rfg.grid[Nx - Nd] - rfg.grid[Nx - Nd + 1] )
    fdAd = fdAd / ( rfg.grid[Nx - Nd + 1] - rfg.grid[Nx - Nd - 1] )
    #@show fdAd

    B = zeros(F64, NCd, NCd)
    Ps = zeros(F64, NCd, Nx)
    Pd = zeros(F64, 4 * Nd, 4 * Nd - 1)

	B[1,1] = 3.0
	B[1,2] = 2.0
	B[1,3] = 1.0
	B[2,1] = 1.0
	B[2,2] = 1.0
	B[2,3] = 1.0
	
	Ps[1,Nx - Nd - 1] = -fdAd
	Ps[1,Nx - Nd + 1] = +fdAd
	Ps[2,Nx - Nd + 0] = +1.0
	Ps[2,Nx - Nd + 1] = -1.0
	
	Pd[1,1] = 1.0
	Pd[2,2] = 1.0
	Pd[3,3] = 1.0
	Pd[4,NCd + 1] = 1.0

    for j = 1:Nd-2
        B[3*j+0,3*j+0] = -RDd[j]
        B[3*j+0,3*j+1] = 3.0
        B[3*j+0,3*j+2] = 2.0
        B[3*j+0,3*j+3] = 1.0
            
        B[3*j+1,3*j-1] = -2 * (RDd[j])^2
        B[3*j+1,3*j+1] = 6.0
        B[3*j+1,3*j+2] = 2.0
            
        B[3*j+2,3*j+1] = 1.0
        B[3*j+2,3*j+2] = 1.0
        B[3*j+2,3*j+3] = 1.0
            
        Ps[3*j+2,Nx-Nd+j+0] = +1.0
        Ps[3*j+2,Nx-Nd+j+1] = -1.0
            
        Pd[4*j+1,3*j+1] = 1.0
        Pd[4*j+2,3*j+2] = 1.0
        Pd[4*j+3,3*j+3] = 1.0
        Pd[4*j+4,NCd+j+1] = 1.0
    end

    j = Nd - 1
	B[3*j+0,3*j+0] = -RDd[j]
	B[3*j+0,3*j+1] = 3.0
	B[3*j+0,3*j+2] = 2.0
	
	B[3*j+1,3*j-1] = -2.0 * (RDd[j])^2
	B[3*j+1,3*j+1] = 6.0
	B[3*j+1,3*j+2] = 2.0
	
	B[3*j+2,3*j+1] = 1.0
	B[3*j+2,3*j+2] = 1.0
	
	Ps[3*j+2,Nx-Nd+j+0] = +1.0
	Ps[3*j+2,Nx-Nd+j+1] = -1.0
	
	Pd[4*j+1,3*j+1] = 1.0
	Pd[4*j+2,3*j+2] = 1.0
	Pd[4*j+4,NCd+j+1] = 1.0

    IB = Matrix{F64}(I, NCd, NCd)
    invB = B \ IB
    #@show invB

    IA = Matrix{F64}(I, Nx, Nx)
    PA = IA[Nx-Nd+1:Nx,1:Nx]
    Ld = vcat(invB * Ps, PA)
    #@show Ld
    Md = Pd * Ld
    #@show Md
    #@show size(Md)

    return Md    
end
