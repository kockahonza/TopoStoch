using DrWatson
@quickactivate "TopoStochSim"

using ComplexAllostery
using ComplexAllostery.Model2_5_C2
using GillespieSim

using GLMakie
using HDF5
using Dates

function m1(N=2; run=nothing)
    ca = makeCAGM(N, 1;
        vars=vars_simplified(; newrs=:ss1, energies=:nocb, concentrations=:ss1)
    )
    @show get_variables(ca)
    # variables should be eb, cATP, cP
    cca = substitute_to_float(ca, [3, 5, 0.0001])

    if !isnothing(run)
        dirpath = datadir("scratch", "simplegil.h5")
        mkpath(dirname(dirpath))

        gs = Gillespie(cca, 0.0, run; rng=100)
        saver = SimpleH5Saver(cca, h5open(dirpath, "w"))

        run_gillespie!(gs;
            saver,
            exit_steps=10000,
            exit_real_time=Second(30)
        )

        close(saver)

        cca, dirpath
    else
        cca
    end
end

function load_from_t1(start_node=1)
    ca = makeCAGM(4, 1;
        vars=vars_simplified(; newrs=:ss1, energies=:nocb, concentrations=:ss1)
    )
    @show get_variables(ca)
    # variables should be eb, cATP, cP
    cca = substitute_to_float(ca, [3, 5, 0.0001])

    filename = projectdir(f"cluster_env/runs/gil_firstrun_localtest/gil1_extreme_3_5_1e-4/start_{start_node}.h5")

    cca, h5open(filename)
end
