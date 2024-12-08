using DrWatson
@quickactivate "cluster_env"

using ComplexAllostery
using ComplexAllostery.Model2_5_C2
using GillespieSim

using Dates

function run_extreme(N)
    ca = makeCAGM(4, 1;
        vars=vars_simplified(; newrs=:ss1, energies=:nocb, concentrations=:ss1)
    )
    @show get_variables(ca)
    cca = substitute_to_float(ca, [3, 5, 0.0001])

    run_full_gillespie_ensemblesim(cca, datadir("gil1_extreme_3_5_1e-4");
        override=true,
        exit_steps=1000000,
        # exit_real_time=Millisecond(100)
        exit_real_time=Hour(14)
    )
end
run_extreme(4)
