using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gm_ca.jl"))

################################################################################
# The simple case, C=2 and B mostly 1
################################################################################
function make_v2_5(N, B; edge_t=:full)
    ca = ComplexAllosteryGM(N, 2, B;
        symmetry=Loop(),
        energy_matrices=make_em_sym(B),
        version=2.5
    )

    if edge_t == :full
        add_edges_full!(ca)
    else
        throw(ArgumentError(f"edge_t \"{edge_t}\" not recognised"))
    end

    ca
end

function add_edge_incweight!(g, u, v, w)
    cur_weight = get_weight(g, u, v)
    add_edge!(g, u, v, cur_weight + w)
end

# define a bunch of symbolic variables so I don't have to keep redoing this
function get_rs()
    r12s = Symbolics.variables(:r, 1:2, 1:2)
    rs = [r12s[i, :] for i in 1:2]
    [rs; Symbolics.variable(:r, 3)]
end
get_thetas() = Symbolics.variables(:θ, 1:3, 1:2)
get_kT() = Symbolics.variable(:kT)
get_chem_energies() = Symbolics.variable.([:εP, :εATP, :εADP])
get_concetrations() = Symbolics.variable.([:cP, :cATP, :cADP])

function add_edges_full!(ca::ComplexAllosteryGM)
    # Declare symbolic variables here
    rs = get_rs()
    thetas = get_thetas()
    kT = get_kT()
    chem_energies = get_chem_energies()
    concetrations = get_concetrations()

    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in (state.conformations[i]+1):ca.C
                new_state = copy(state)
                new_state.conformations[i] = new_c
                new_vertex = statetoi(new_state, ca)

                e_i_state = calc_monomer_energy(state, i, ca) + calc_interaction_energy(state, i, ca)
                e_i_new_state = calc_monomer_energy(new_state, i, ca) + calc_interaction_energy(new_state, i, ca)

                # forward rate (where c_i increases)
                exp_term_f = thetas[3, 1] * e_i_state - (1 - thetas[3, 2]) * e_i_new_state
                rf = rs[3] * exp(exp_term_f / kT)
                add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                # backward rate (where c_i decreases)
                exp_term_b = thetas[3, 2] * e_i_new_state - (1 - thetas[3, 1]) * e_i_state
                rb = rs[3] * exp(exp_term_b / kT)
                add_edge_incweight!(ca.graph, new_vertex, vertex, rb)
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                new_vertex = statetoi(new_state, ca)

                e_i_state = calc_monomer_energy(state, i, ca)
                e_i_new_state = calc_monomer_energy(new_state, i, ca)

                i_conf = state.conformations[i]

                # P + S <-> NS reaction
                begin
                    # forward rate
                    exp_term_f = thetas[1, 1] * (e_i_state + chem_energies[1]) - (1 - thetas[1, 2]) * e_i_new_state
                    rf = rs[1][i_conf] * concetrations[1] * exp(exp_term_f / kT)
                    add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[1, 2] * e_i_new_state - (1 - thetas[1, 1]) * (e_i_state + chem_energies[1])
                    rb = rs[1][i_conf] * exp(exp_term_b / kT)
                    add_edge_incweight!(ca.graph, new_vertex, vertex, rb)
                end

                # ATP + S <-> ADP + NS reaction
                begin
                    # forward rate
                    exp_term_f = thetas[2, 1] * (e_i_state + chem_energies[2]) - (1 - thetas[2, 2]) * (e_i_new_state + chem_energies[3])
                    rf = rs[2][i_conf] * concetrations[2] * exp(exp_term_f / kT)
                    add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[2, 2] * (e_i_new_state + chem_energies[3]) - (1 - thetas[2, 1]) * (e_i_state + chem_energies[1])
                    rb = rs[2][i_conf] * concetrations[3] * exp(exp_term_b / kT)
                    add_edge_incweight!(ca.graph, new_vertex, vertex, rb)
                end
            end
        end
    end
end

function fsimp(ca::ComplexAllosteryGM; which=:both)
    if which == :both
        which = [:thetas, :chem_energies]
    elseif isa(which, Symbol)
        which = [which]
    end

    terms = Dict{Num,Num}()
    if :thetas in which
        for theta in get_thetas()
            terms[theta] = 1.0
        end
    end
    if :chem_energies in which
        for chem_energy in get_chem_energies()
            terms[chem_energy] = 0.0
        end
    end

    terms[get_kT()] = 1.0

    ssubstitute(ca, terms)
end

function frsimp(ca::ComplexAllosteryGM; do_fsimp=true, kwargs...)
    rs = get_rs()
    terms = Dict{Num,Num}()

    terms[rs[1][1]] = 1.0
    terms[rs[1][2]] = 1.0

    terms[rs[2][1]] = 1.0
    terms[rs[2][2]] = 0.0

    terms[rs[3]] = 1.0

    if do_fsimp
        ca = fsimp(ca; kwargs...)
    end
    ssubstitute(ca, terms)
end

function frsimp2(ca::ComplexAllosteryGM; do_fsimp=true, kwargs...)
    rs = get_rs()
    newrs = Symbolics.variables(:r, 1:3)
    alpha = Symbolics.variable(:α)

    terms = Dict{Num,Num}()

    terms[rs[1][1]] = newrs[1]
    terms[rs[1][2]] = newrs[1]

    terms[rs[2][1]] = newrs[2]
    terms[rs[2][2]] = newrs[2] * alpha

    terms[rs[3]] = newrs[3]

    if do_fsimp
        ca = fsimp(ca; kwargs...)
    end
    ssubstitute(ca, terms), vcat(newrs, alpha)
end

function int_plot_1(ca::ComplexAllosteryGM, args...; init_scen=:eonly, kwargs...)
    if ca.version != 2.5
        throw(ArgumentError("this method only works for version 2.5"))
    end
    sca, frterms = frsimp2(ca)
    variables = vcat(get_em_vars(), get_concetrations(), frterms)
    ranges = [
        -1.0:0.1:10.0,
        0.0:0.1:5.0,
        0.0:0.1:10.0,
        0.0:0.1:5.0,
        0.0:0.1:5.0,
        0.0:0.01:2.0,
        0.0:0.01:1.0,
        0.0:0.1:5.0,
        0.0:0.05:5.0,
        0.0:0.01:1.5
    ]
    if init_scen == :eonly
        startvalues = [3.0, 0.5, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    elseif init_scen == :int1
        startvalues = [3.0, 0.5, 0.2, 1.0, 5.0, 0.05, 0.25, 2.0, 0.25, 0.0]
    end
    int_plot_ca(sca, args...; variables, ranges, startvalues, kwargs...)
end


function get_eq_subs1()
    terms = Dict{Num,Num}()
    es = get_chem_energies()
    terms[es[1]] = es[2] - es[3]
    cs = get_concetrations()
    terms[cs[1]] = cs[2] / cs[3]
    terms
end
