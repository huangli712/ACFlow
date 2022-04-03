#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

welcome()
overview()
read_param()
ω, 𝐺 = read_data!()
Aom = som_run(ω, 𝐺)
som_output(Aom)
