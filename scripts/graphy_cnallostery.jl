using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("GraphModels.jl"))

using Graphs, MetaGraphsNext, SimpleWeightedGraphs
using GLMakie, GraphMakie, NetworkLayout
using StatsBase
using Symbolics
using StaticArrays, DataFrames
using PrettyTables
using PyFormattedStrings

import Base: copy, broadcastable, show

################################################################################
# Needed combinatoric utilities
################################################################################

"""
Returns true if `perm` has length `N` and all its elements are between `base`(inclusive) and `base` + `E`(exclusive).
"""
function validperm(perm::AbstractVector{<:Integer}, E, N; base=1)
    if length(perm) != N
        return false
    end
    for e in perm
        if (e < base) || (e >= base + E)
            return false
        end
    end
    return true
end

"""
`perm` is assumed to be a valid permutation and `E` is as in `validperm`.
Returns in the index of the given permutation.
"""
function permtoi(perm::AbstractVector{<:Integer}, E; base=1)
    1 + sum((ei - base) * E^(i - 1) for (i, ei) in enumerate(perm))
end

"""
Returns the `i`th permutation of length `N` with each element being an
integer between 1 and `E`. The ordering is the same order as in `permtoi`.
"""
function itoperm(i, E, N; base=1)
    perm = Vector{Int}(undef, N)
    for j in N:-1:1
        factor = E^(j - 1)
        perm[j] = div(i - 1, factor)
        i -= perm[j] * factor
        perm[j] += base
    end
    perm
end

################################################################################
# Complex Allostery types
################################################################################
struct CAState{T}
    conformations::Vector{T}
    occupations::Vector{T}
end
function copy((; conformations, occupations)::CAState)
    CAState(copy(conformations), copy(occupations))
end
show(io::IO, st::CAState) = print(io, f"({st.conformations}, {st.occupations})")

abstract type Symmetry end
struct Chain <: Symmetry end
struct Loop <: Symmetry end
broadcastable(s::Symmetry) = Ref(s)

struct EnergyMatrices{F<:Number}
    monomer::Matrix{F}
    interactions::Matrix{F}
end
broadcastable(em::EnergyMatrices) = Ref(em)

mutable struct ComplexAllosteryGM{S<:Symmetry,F<:Number,T<:Integer} <: AbstractGraphModel
    N::T # Number of monomers
    C::T # Number of possible conformational states per monomer
    B::T # Number of ligand binding sites per monomer

    energy_matrices::EnergyMatrices

    graph::SimpleWeightedDiGraph{T,F}
    function ComplexAllosteryGM(N::T, C::T, B::T; _::S=Chain(), energy_matrices=nothing, sym_graph=false) where {T<:Integer,S}
        Ftype = sym_graph ? Num : Float64

        graph = SimpleWeightedDiGraph{T,Ftype}(numstates(N, C, B))

        if isnothing(energy_matrices)
            energy_matrices = make_EM_sym_general(C, B)
        else
            if !validEM(energy_matrices, C, B)
                throw(ArgumentError(f"the provided energy matrices do not have the right shape for the given system"))
            end
        end
        new{S,Ftype,T}(N, C, B, energy_matrices, graph)
    end
end
broadcastable(ca::ComplexAllosteryGM) = Ref(ca)
graph(ca::ComplexAllosteryGM) = ca.graph

# Utils state related functions
function numstates(N, C, B)
    (C * (B + 1))^N
end
numstates(ca::ComplexAllosteryGM) = numstates(ca.N, ca.C, ca.B)

function validstate(st::CAState, N, C, B)
    validperm(st.conformations, C, N) && validperm(st.occupations, B + 1, N; base=0)
end
validstate(st::CAState, ca::ComplexAllosteryGM) = validstate(st, ca.N, ca.C, ca.B)

function statetoi(st::CAState, N, C, B)
    permtoi(st.conformations, C) + (permtoi(st.occupations, B + 1; base=0) - 1) * (C^N)
end
statetoi(st::CAState, ca::ComplexAllosteryGM) = statetoi(st, ca.N, ca.C, ca.B)

function itostate(i, N, C, B)
    (iocc, iconf) = divrem(i - 1, C^N)
    CAState(itoperm(iconf + 1, C, N), itoperm(iocc + 1, B + 1, N; base=0))
end
itostate(i, ca::ComplexAllosteryGM) = itostate(i, ca.N, ca.C, ca.B)

allstates(N, C, B) = (itostate(i, N, C, B) for i in 1:numstates(N, C, B))
allstates(ca::ComplexAllosteryGM) = allstates(ca.N, ca.C, ca.B)

################################################################################
# Dealing with energy and equilibrium probability calculations
################################################################################
# Setting up the energy matrices
function make_EM_sym_general(C, B)
    monomer = Matrix{Num}(undef, C, B + 1)
    monomer[1, 1] = 0
    monomer[2:C, 1] = Symbolics.variables(:εr, 1:C-1)
    monomer[:, 2:B+1] = Symbolics.variables(:εl, 1:C, 1:B)

    interactions = Matrix{Num}(undef, C, C)
    terms = Symbolics.variables(:εb, 1:div((C - 1) * C, 2))
    term_i = 1
    for i in 1:C
        interactions[i, i] = 0
        for j in i+1:C
            interactions[i, j] = terms[term_i]
            interactions[j, i] = terms[term_i]
            term_i += 1
        end
    end

    EnergyMatrices(monomer, interactions)
end

function make_EM_sym_C2(B)
    monomer = Matrix{Num}(undef, 2, B + 1)
    monomer[1, 1] = 0
    monomer[1, 2:B+1] = Symbolics.variables(:εt, 1:B)
    monomer[2, 1:B+1] = Symbolics.variables(:εr, 0:B)

    eb = Symbolics.variable(:εb)
    interactions = [0 eb; eb 0]
    EnergyMatrices(monomer, interactions)
end

function make_EM_sym_C2_simpler(B)
    monomer = Matrix{Num}(undef, 2, B + 1)
    @variables ε_t, Δε_r
    monomer[1, 1] = 0
    monomer[2, 1] = Δε_r
    monomer[1, 2:B+1] .= ε_t .* collect(1:B)
    monomer[2, 2:B+1] .= monomer[1, 2:B+1] .- Δε_r

    eb = Symbolics.variable(:ε_b)
    interactions = [0 eb; eb 0]
    EnergyMatrices(monomer, interactions)
end

function validEM(em::EnergyMatrices, C, B)
    (size(em.monomer) == (C, B + 1)) && (size(em.interactions) == (C, C))
end
validEM(em::EnergyMatrices, ca::ComplexAllosteryGM) = validEM(em, ca.C, ca.B)

function get_variables(em::EnergyMatrices)
    unique(Iterators.flatten(Symbolics.get_variables.(Iterators.flatten((em.monomer, em.interactions)))))
end

# Simple calculator functions
function calc_energy(st::CAState, em::EnergyMatrices, ::Chain)
    energy = 0
    for i in 1:length(st.conformations)
        energy += em.monomer[st.conformations[i], st.occupations[i]+1]
    end
    energy += 0.5 * em.interactions[st.conformations[1], st.conformations[2]]
    for i in 2:length(st.conformations)-1
        energy += 0.5 * (
            em.interactions[st.conformations[i-1], st.conformations[i]] +
            em.interactions[st.conformations[i], st.conformations[i+1]]
        )
    end
    energy += 0.5 * em.interactions[st.conformations[end-1], st.conformations[end]]
    energy
end
function calc_energy(st::CAState, em::EnergyMatrices, ::Loop)
    energy = 0
    for i in 1:length(st.conformations)
        energy += em.monomer[st.conformations[i], st.occupations[i]+1]
    end
    energy += 0.5 * (
        em.interactions[st.conformations[end], st.conformations[1]] +
        em.interactions[st.conformations[1], st.conformations[2]]
    )
    for i in 2:length(st.conformations)-1
        energy += 0.5 * (
            em.interactions[st.conformations[i-1], st.conformations[i]] +
            em.interactions[st.conformations[i], st.conformations[i+1]]
        )
    end
    energy += 0.5 * (
        em.interactions[st.conformations[end-1], st.conformations[end]] +
        em.interactions[st.conformations[end-1], st.conformations[end]]
    )
    energy
end
function calc_energy(st::CAState, ca::ComplexAllosteryGM{S}) where {S}
    calc_energy(st, ca.energy_matrices, S())
end
calc_numligands(st::CAState) = sum(st.occupations)
function calc_gibbs_factor(st::CAState, em::EnergyMatrices, symmetry)
    @variables kT, μ
    exp(-(calc_energy(st, em, symmetry) - μ * calc_numligands(st)) / kT)
end
calc_gibbs_factor(st::CAState, ca::ComplexAllosteryGM{S}) where {S} = calc_gibbs_factor(st, ca.energy_matrices, S())

# Calculating equilibrium observables
function calc_all_gibbs_factors(ca::ComplexAllosteryGM{S}) where {S}
    calc_gibbs_factor.(allstates(ca), ca)
end
function calc_partition_function(ca::ComplexAllosteryGM)
    sum(calc_all_gibbs_factors(ca))
end
function calc_avg_numligands(ca::ComplexAllosteryGM)
    factors = calc_all_gibbs_factors(ca)
    sum(calc_numligands.(allstates(ca)) .* factors) / sum(factors)
end
function calc_avg_energy(ca::ComplexAllosteryGM)
    factors = calc_all_gibbs_factors(ca)
    sum(calc_energy.(allstates(ca), ca) .* factors) / sum(factors)
end

################################################################################
# Dealing with transitions
################################################################################
# Adding some basic links without meaningful weights
function add_edges_base!(ca::ComplexAllosteryGM)
    # Just use a constant for all the rates,disregarding detailed balance
    universal_rate = 0.75 / (ca.C + 1) # This means that the chance of not doing anything is ~0.25
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
            end
        end
    end
end

function add_edges_sym!(ca::ComplexAllosteryGM)
    rates_i = 1
    function add_new_rate(v1, v2)
        if (v1 < v2)
            rate = Symbolics.variable(:r, rates_i)
            add_edge!(ca.graph, v1, v2, rate)
            add_edge!(ca.graph, v2, v1, rate)
            rates_i += 1
        end
    end
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_new_rate(vertex, statetoi(new_state, ca))
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
        end
    end
end

################################################################################
# The simple case, C=2 and B mostly 1
################################################################################

function make_simple(N, B)
    ca = ComplexAllosteryGM(N, 2, B; energy_matrices=make_EM_sym_C2_simpler(B), sym_graph=true)

    # Label occupancy changes using r_(ci)+-
    occ_increase_rates = Symbolics.variables(:a, 1:2)
    occ_decrease_rates = Symbolics.variables(:m, 1:2)
    # Label conformational changes via increasing indices or ri
    rates_i = 1
    function add_new_indexed_rate(v1, v2)
        if (v1 < v2)
            rate = Symbolics.variable(:ri, rates_i)
            add_edge!(ca.graph, v1, v2, rate)
            add_edge!(ca.graph, v2, v1, rate)
            rates_i += 1
        end
    end
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_new_indexed_rate(vertex, statetoi(new_state, ca))
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), occ_decrease_rates[state.conformations[i]])
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), occ_increase_rates[state.conformations[i]])
            end
        end
    end

    ca
end

function find_a1(ca::ComplexAllosteryGM)
    el = Symbolics.variable(:a, 1)
    am = adjacency_matrix(ca.graph)

    cis = findall(x -> isequal(x, el), am)

    equal_items = []
    # push!(equal_items, el / Symbolics.variable(:m, 1))
    for ci in cis
        (i, j) = Tuple(ci)
        push!(equal_items, calc_gibbs_factor(itostate(j, ca), ca) / calc_gibbs_factor(itostate(i, ca), ca))
    end

    equal_items
end

################################################################################
# Util and run function
################################################################################
# Plotting, also used by simGM
function plotGM(ca::ComplexAllosteryGM{S,F}, args...;
    interactive=false, sub_layout=Shell(), z_scale=1.0, layout=nothing, kwargs...
) where {S,F}
    @inline function get_adj_matrix()
        if F == Num
            map(x -> if iszero(x)
                    0.0
                else
                    1.0
                end, adjacency_matrix(ca.graph))
        else
            adjacency_matrix(ca.graph)
        end
    end

    if isnothing(layout)
        layout = [(x, y, z * z_scale) for ((x, y), z) in zip(sub_layout(get_adj_matrix()), calc_numligands.(allstates(ca)))]
    elseif isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(get_adj_matrix())
    end

    node_color = [:snow3 for _ in 1:numstates(ca)]
    node_size = [10 for _ in 1:numstates(ca)]

    fap = graphplot(ca.graph, args...; layout, node_color, node_size, curve_distance_usage=false, kwargs...)

    if interactive
        make_interactive(fap)
    end

    fap
end

function plot_CA_sym(ca::ComplexAllosteryGM, args...; kwargs...)
    plotGM(ca, args...; nlabels=repr.(allstates(ca)), elabels=repr.(weight.(edges(ca.graph))), elabels_shift=0.3, kwargs...)
end

function eq_stats_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    variables = Num.(get_variables(ca.energy_matrices))
    @variables μ, kT
    append!(variables, μ)

    sliders = [(label=repr(var), range=0.0:0.01:10.0, startvalue=1.0) for var in variables]
    sg = SliderGrid(fig[1, 1], sliders...)
    variable_observables = [x.value for x in sg.sliders]

    overall_stats = GridLayout()
    fig[1, 2] = overall_stats
    colsize!(fig.layout, 2, 30)

    partition_function_f = build_function(substitute(calc_partition_function(ca), (kT => 1)), variables; expression=Val(false))

    avg_numligands_f = build_function(substitute(calc_avg_numligands(ca), (kT => 1)), variables; expression=Val(false))
    avg_numligands_o = lift((args...) -> f"N={avg_numligands_f(args):.3g}", variable_observables...)
    Label(overall_stats[1, 1], avg_numligands_o)

    avg_energy_f = build_function(substitute(calc_avg_energy(ca), (kT => 1)), variables; expression=Val(false))
    avg_energy_o = lift((args...) -> f"E={avg_energy_f(args):.3g}", variable_observables...)
    Label(overall_stats[2, 1], avg_energy_o)

    # Do the barplot
    state_reordering = sortperm(calc_numligands.(allstates(ca)); alg=Base.Sort.DEFAULT_STABLE)
    state_labels = repr.(allstates(ca))[state_reordering]
    ax = Axis(fig[2, :];
        xticks=(1:numstates(ca), state_labels),
        xticklabelrotation=pi / 5
    )
    ylims!(ax, 0.0, 0.5)

    gfs = substitute.(calc_gibbs_factor.(allstates(ca), ca), (kT => 1.0))
    gfs_fs = build_function.(gfs, Ref(variables); expression=Val(false))

    probabilities = lift(variable_observables...) do (vars...)
        pfval = partition_function_f(vars)
        [gfs_fs[i](vars) / pfval for i in 1:numstates(ca)][state_reordering]
    end

    bp = barplot!(ax, 1:numstates(ca), probabilities)

    Makie.FigureAxisPlot(fig, ax, bp)
end

function eq_coopbinding_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    energy_vars = get_variables(ca.energy_matrices)
    @variables μ, εsol, kT, c

    slider_vars = [energy_vars; εsol]

    sliders = [(label=repr(var), range=0.0:0.01:10.0, startvalue=1.0) for var in slider_vars]
    sg = SliderGrid(fig[1, 1], sliders...)
    slider_observables = [x.value for x in sg.sliders]

    ax = Axis(fig[2, :])

    prepped_expression = substitute(calc_avg_numligands(ca), Dict(
        kT => 1,
        μ => εsol + log(c)
    ))
    avg_numligands_f = build_function(prepped_expression, [energy_vars; εsol; c]; expression=Val(false))

    conc_range = LinRange(0, 10, 300)

    avg_numligands = lift(slider_observables...) do (vars...)
        [avg_numligands_f((vars..., cval)) for cval in conc_range]
    end

    sc = scatter!(ax, conc_range, avg_numligands)

    Makie.FigureAxisPlot(fig, ax, sc)
end

function save_adj_matrix(ca::ComplexAllosteryGM; name=savename("adjmat", ca), short=false)
    file = open(plotsdir(savename("adjmat", ca, "table")), "w")

    row_names = repr.(1:numstates(ca)) .* " / " .* repr.(allstates(ca))
    modified_matrix = [row_names Matrix{Any}(adjacency_matrix(ca.graph))]
    header = if short
        [""; repr.(1:numstates(ca))]
    else
        [""; repr.(allstates(ca))]
    end

    pretty_table(file, modified_matrix; header)
    close(file)
end

function save_gibbs_factors(ca::ComplexAllosteryGM; name=savename("adjmat", ca))
    file = open(plotsdir(savename("gibbsfs", ca, "table")), "w")
    data = [repr.(allstates(ca)) repr.(calc_gibbs_factor.(allstates(ca), ca))]
    pretty_table(file, data; header=["state", "gibbs factor"])
    close(file)
end
