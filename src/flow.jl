function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
end

function read_data()
    finput = get_c("finput")
    ngrid = get_c("ngrid")
    grid = get_c("grid")
    @show finput, ngrid, grid

    if grid == "matsubara"
        return read_freq_data(finput, ngrid)
    else
        return read_time_data(finput, ngrid)
    end
end

function solve()
    solver = get_c("solver")
    
    @cswitch solver begin
        @case "MaxEnt"
            break

        @case "StochOM"
            break

        @case "StochAC"
            break

        @default
            sorry()
            break
    end
end
