using DrWatson
@quickactivate "cluster_env"

using ComplexAllostery
using ComplexAllostery.Model2_5_C2
using GillespieSim

using Dates

function run_extreme(N=4; test=false)
    ca = makeCAGM(N, 1;
        vars=vars_simplified(; newrs=:ss1, energies=:nocb, concentrations=:ss1)
    )
    @show get_variables(ca)
    # variables should be eb, cATP, cP
    cca = substitute_to_float(ca, [3, 5, 0.0001])

    if test
        exit_real_time = Millisecond(500)
    else
        exit_real_time = Hour(6)
    end

    run_full_gillespie_ensemblesim(cca, "./gil1_extreme_3_5_1e-4";
        override=true,
        exit_steps=100000,
        exit_real_time,
        rng=100
    )
end
# run_extreme(4)
