using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("gm_complexallostery.jl"))

################################################################################
# The simple case, C=2 and B mostly 1
################################################################################
function make_simple(N, B; edge_t=:EL1)
    ca = ComplexAllosteryGM(N, 2, B; energy_matrices=make_EM_sym_C2(B))

    if edge_t == :EL1
        add_edges_simple_EL1!(ca)
    elseif edge_t == :EL2
        add_edges_simple_EL2!(ca)
    else
        throw(ArgumentError(f"edge_t \"{edge_t}\" not recognised"))
    end

    ca
end

"""
Energy landscape inspired symbolic edges that use energy height barriers
and and some physicsy arguments done on pen and paper. Maybe works?
"""
function add_edges_simple_EL1!(ca::ComplexAllosteryGM)
    rate_occ_change = Symbolics.variable(:r_B)
    rate_conf_change = Symbolics.variable(:r_C)
    @variables ε_t, Δε_r, ε_b, μ, kT
    fact_1_dec = exp(-(-ε_t + μ) / kT)
    fact_2_inc_0 = exp(-(-Δε_r) / kT)
    fact_2_inc = exp(-(Δε_r) / kT)
    fact_2_dec = exp(-(-ε_t + μ + Δε_r) / kT)

    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    exponent_term = Num(0)
                    if state.occupations[i] != 0
                        exponent_term += ε_t
                        if state.conformations[i] == 2
                            exponent_term -= Δε_r
                        end
                    else
                        if state.conformations[i] == 2
                            exponent_term += Δε_r
                        end
                    end
                    exponent_term += calc_neighbors_energy(state, i, ca)
                    rate = rate_conf_change / exp(-exponent_term / kT)
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate
                    )
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            cur_con = state.conformations[i]
            if cur_con == 1
                if cur_occ < ca.B
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change
                    )
                end
                if cur_occ > 0
                    new_state = copy(state)
                    new_state.occupations[i] -= 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_1_dec
                    )
                end
            elseif cur_con == 2
                if cur_occ == 0
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_inc_0
                    )
                elseif cur_occ < ca.B
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_inc
                    )
                end
                if cur_occ > 0
                    new_state = copy(state)
                    new_state.occupations[i] -= 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_dec
                    )
                end
            else
                throw(ErrorException("kaka"))
            end
        end
    end
end

"""
Energy landscape inspired symbolic edges that use energy height barriers
and the actual gibbs factors for getting the rates. This however could very
well lead to rates > 1. Hence not ready!
"""
function add_edges_simple_EL2!(ca::ComplexAllosteryGM)
    gfs = calc_gibbs_factor.(allstates(ca), ca)

    # Label occupancy changes using r_(ci)+-
    rate_occ_change = Symbolics.variable(:r)
    @variables ε_t, Δε_r, ε_b, μ, kT

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
                    # add_new_indexed_rate(vertex, statetoi(new_state, ca))
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                rate = rate_occ_change * exp(-cur_occ * (ε_t - μ) / kT) / gfs[vertex]
                add_edge!(
                    ca.graph,
                    vertex,
                    statetoi(new_state, ca),
                    simplify(rate; rewriter=get_rewriter())
                )
            end
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                rate = rate_occ_change * exp(-(cur_occ - 1) * (ε_t - μ) / kT) / gfs[vertex]
                add_edge!(
                    ca.graph,
                    vertex,
                    statetoi(new_state, ca),
                    simplify(rate; rewriter=get_rewriter())
                )
            end
        end
    end
end

function terms_simple(mu, et, der, eb, rB, rC, kT_=1.0)
    @variables μ, ε_t, Δε_r, ε_b, r_B, r_C, kT
    Dict(
        μ => mu,
        ε_t => et,
        Δε_r => der,
        ε_b => eb,
        r_B => rB,
        r_C => rC,
        kT => kT_
    )
end
terms_simple0() = terms_simple(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

function terms_simple1(et, eb_option, mu_option)
    eb = if eb_option == 1
        0.0
    elseif eb_option == 2
        0.5 * et
    else
        throw(ArgumentError("eb_option must be 1 or 2"))
    end
    mu = if mu_option == 1
        0.3 * et
    elseif mu_option == 2
        0.8 * et
    elseif mu_option == 3
        1.2 * et
    else
        throw(ArgumentError("mu_option must be 1, 2 or 3"))
    end
    terms_simple(mu, et, 0.2 * et, eb, exp(-et), exp(-2.5 * et), 1.0)
end

function run_sample_casimple()
    ca = substitute_to_float(make_simple(4, 1), terms_simple1(4.0, 2, 3))
    simGM(ca, 2000; num_vertices=10, layout=(:NRbs3,), delay=0.05)
end
