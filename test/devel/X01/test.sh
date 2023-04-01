cp Aout.data Aout.data.0

python test.py > fmesh.inp
../../../util/acrun.jl ac.toml
cp Aout.data Aout.data.1

python test.py > fmesh.inp
../../../util/acrun.jl ac.toml
cp Aout.data Aout.data.2

python test.py > fmesh.inp
../../../util/acrun.jl ac.toml
cp Aout.data Aout.data.3

python test.py > fmesh.inp
../../../util/acrun.jl ac.toml
cp Aout.data Aout.data.4

python test.py > fmesh.inp
../../../util/acrun.jl ac.toml
cp Aout.data Aout.data.5
