#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/01/30
#

function make_model(m::AbstractMesh)
    model = get_c("model")
    if model == "flat"
        return make_flat_model(m)
    else
        return make_gaussian_model(m)
    end
end

function make_flat_model(um::UniformMesh)
    model = ones(F64, um.nmesh)
    norm = dot(um.weight, model)
    model = model ./ norm
    return model
end

function make_flat_model(num::NonUniformMesh)
end

function make_gaussian_model(um::UniformMesh)
end

function make_gaussian_model(num::NonUniformMesh)
    model = exp.(-(num.mesh / 20.0) .^ 6.0)
    norm = dot(num.weight, model)
    model = model ./ norm
    return model
end