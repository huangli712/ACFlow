for i in {1..100}
do
    echo $i
    ../../util/acrun.jl ac.toml
done
