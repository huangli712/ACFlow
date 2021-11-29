#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/29
#

function calc_mesh(𝑀::MomentsData)
    # Ratio of main spectral range and standard deviation of spectrum
    R_SW_STD_OMEGA = 3.0
    
    # Minimum ratio of standard deviation and frequency step
    RMIN_SW_DW = 50.0

    # Total real frequency range
    F_W_RANGE = 20.0

    SC = 𝑀.𝑀₁ / 𝑀.𝑀₀
    SW = sqrt(𝑀.𝑀₂ / 𝑀.𝑀₀ - SC^2) * R_SW_STD_OMEGA
    wl = SC - SW / 2
    wr = SC + SW / 2
    origin = SC
    dw = SW / R_SW_STD_OMEGA / RMIN_SW_DW

    println("SC: ", SC)
    println("SW: ", SW)
    println("wl: ", wl)
    println("wr: ", wr)
    println("origin: ", origin)
    println("dw: ", dw)

    nwr = round(I64, (wr - origin) / dw)
    wr = origin + nwr * dw
    wcr = collect(range(origin, wr, length =  nwr + 1))

    nwl = round(I64, (origin - wl) / dw)
    wl = origin - nwl * dw
    wcl = collect(range(wl, origin, length = nwl + 1)) 
    #@show nwl, nwr, wl, wr
    #@show wcr
    #@show wcl

    wc = [wcl[1:end-1]; wcr]
    #@show length(wc), wc

    wmin = SC - F_W_RANGE * SW / 2.0
    w0l = wl + sqrt(dw * (wl - wmin))
    #@show w0l
    nul = ceil((w0l - wl) / dw)
    w0l = wl + nul * dw
    dul = dw / (wl - dw - w0l) / (wl - w0l)
    wmin = -1.0 / dul + w0l
    #@show wmin, w0l, nul, dul
    ul = F64[]
    for i = 1:nul
        push!(ul, -dul * i)
    end
    #for i = 1:length(ul)
    #    @show i, ul[i]
    #end

    wmax = SC + F_W_RANGE * SW / 2.0
    w0r = wr - sqrt(dw * (wmax - wr))
    #@show w0r
    nur = ceil((wr - w0r) / dw)
    w0r = wr - nur * dw
    dur = dw / (wr + dw - w0r) / (wr - w0r)
    wmax = 1.0 / dur + w0r
    #@show wmax, w0r, nur, dur
    ur = F64[]
    for i = 1:nur
        insert!(ur, 1, dur * i)
    end
    #for i = 1:length(ur)
    #    @show i, ur[i]
    #end

    w = [1.0 ./ ul .+ w0l; wc; 1.0 ./ ur .+ w0r]
    #@show length(w)
    #for i = 1:length(w)
    #    @show i, w[i]
    #end

    return RealFrequencyGrid(w)
end
