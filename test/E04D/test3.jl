#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using DelimitedFiles
using ACFlow

nfine = 1000000
wmin  = -5.0
wmax  = +5.0
dlm = readdlm("Aout.data")

@show size(dlm)

mesh = LinearMesh(dlm[:,1])
Aout = dlm[:,2]
fmesh = LinearMesh(nfine, wmin, wmax)

allow = zeros(I64, nfine)
#@show allow
for i = 1:nfine
    p = nearest(mesh, i / nfine)
    if Aout[p] > 0.02
        allow[i] = 1
    end
end

@show count(==(1), allow)