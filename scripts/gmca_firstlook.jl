using DrWatson
@quickactivate "TopoStochSim"

using Revise

using GLMakie

using ComplexAllostery
using ComplexAllostery.Model2_5_C2

################################################################################
# 24-12-03 - looking at the cycles in N~3-6 and the pertrubation states
function m1()
    ca = makeCAGM(4, 1; vars=vars_simplified(; newrs=:ss1, energies=false, concentrations=:ss1))
    @show get_variables(ca)
    cca = substitute_to_float(ca, [0, 0, 2, 10, 0.001])

    sc1 = GLMakie.Screen()
    display(sc1, plotcgm(cca; layout=:NTc))

    es = fsteadystates(cca)[2]
    @show findall(isnothing, degenerate_map(es))

    display(plot_ps(cca, evec(es, 1); layout=:NTc))

    ca, cca, es, sc1
end
# ca, cca, es, sc1 = m1()

# 24-12-12 - Just a helper function to quickly make a ca
function m2()
    ca = makeCAGM(4, 1; vars=vars_simplified(; newrs=:ss1, energies=false, concentrations=:ss1))
    @show get_variables(ca)
    cca = substitute_to_float(ca, [0, 0, 2, 10, 0.001])

    ca, cca
end
# ca, cca = m2()
