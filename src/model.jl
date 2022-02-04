#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

function make_flat_model(um::UniformMesh)
    model = ones(F64, um.nmesh)
    norm = dot(um.weight, model)
    model = model ./ norm
    return model
end

function make_flat_model(num::NonUniformMesh)
end

function make_gaussian_model(um::UniformMesh)
    # For test8
    #println("here")
    model = exp.(-0.5 * (um.mesh) .^ 2.0)
    norm = dot(um.weight, model)
    model = model ./ norm
    return model
end

function make_gaussian_model(num::NonUniformMesh)
    # For test5
    #model = exp.(-(num.mesh / 20.0) .^ 6.0)
    #norm = dot(num.weight, model)
    #model = model ./ norm .* 4.77

    # For test6
    model = exp.(-(num.mesh / 8.0) .^ 6.0)
    norm = dot(num.weight, model)
    model = model ./ norm
    return model
end

function make_file_model(am::AbstractMesh)
    model = zeros(F64, am.nmesh)

    open("model.data", "r") do fin
        for i = 1:am.nmesh
            model[i] = parse(F64, line_to_array(fin)[2])
        end
    end

    return model
end