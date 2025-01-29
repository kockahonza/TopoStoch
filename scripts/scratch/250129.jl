using DrWatson
@quickactivate "TopoStochSim"

using Revise

using GLMakie

using ComplexAllostery
using ComplexAllostery.Model3
using GillespieSim

################################################################################
# 25-01-29 - first look at model 3
function m1()
    ca = makeCAGM(4, vars_simplified())
    @show get_variables(ca)

    # Display a concrete sample in a persistent window
    display(GLMakie.Screen(), plotcgm(substitute_to_float(ca, [1, 1, 10, 1, 1, 0.01, 0, 0]); layout=:NTc))

    ca
end
