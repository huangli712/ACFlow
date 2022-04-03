#!/usr/bin/env julia -p 8
  
@everywhere push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

@everywhere using ACFlow

welcome()
overview()
read_param()
ω, 𝐺 = read_data!()
Aom = sum( pmap((x) -> som_run(ω, 𝐺), 1:nworkers()) ) ./ nworkers()
som_output(Aom)
