#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/01/29
#

function solve(rd::RawData)
    solver = get_c("solver")
    
    @cswitch solver begin
        @case "MaxEnt"
            MaxEnt.solve(rd)
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
