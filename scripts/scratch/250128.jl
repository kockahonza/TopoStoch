using DrWatson
@quickactivate "TopoStochSim"

using Revise

using GLMakie

using ComplexAllostery
using ComplexAllostery.Model2_6_C2
using GillespieSim

################################################################################
# 25-01-28 - first look at model 2.6
function m1()
    ca = makeCAGM(4, vars_simplified())
    @show get_variables(ca)

    # Display a concrete sample in a persistent window
    display(GLMakie.Screen(), plotcgm(substitute_to_float(ca, [1, 1, 10, 1, 1, 0.01, 0, 0]); layout=:NTc))

    ca
end

function model2_6_examples_from_meeting()
    ca = makeCAGM(4, vars_simplified())
    # Should be Num[ε_b, r₃, cATP, r₂, r₁, cP, α, cADP]
    @show get_variables(ca)

    # Two disjoint cycles
    display(GLMakie.Screen(), plotcgm(substitute_to_float(ca, [5, 1, 6, 1, 1, 0.0001, 0, 0]); layout=:NTc))
    # One global cycle
    display(GLMakie.Screen(), plotca_macro(substitute_to_float(ca, [5, 10000, 6, 1, 0.1, 0.0001, 0, 0])))
    # Shrinking rectangle
    display(GLMakie.Screen(), plotca_macro(substitute_to_float(ca, [5, 10000, 6, 0.05, 0.1, 0.0001, 0, 0])))
    # A decent triangle
    display(GLMakie.Screen(), plotca_macro(substitute_to_float(ca, [5, 10000, 6, 20, 0.1, 0.0001, 0, 0])))

    # kinda triangle and also multiple overlaid futile cycles
    display(GLMakie.Screen(), plotca_macro(substitute_to_float(ca, [5, 1000, 6, 5, 0.01, 0.0001, 0, 0])))
end
