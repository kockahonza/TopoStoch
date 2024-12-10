using DrWatson
@quickactivate "TopoStochSim"

using Revise

using Ones
using GillespieSim

using GLMakie
using Dates

function m1()
    ogm = make_single_cycle(5, get_c2d_int(1))
    l = Gillespie(ogm, 0.0, 1; rng=100)

    ogm, l
end
# ogm, l = m1()
