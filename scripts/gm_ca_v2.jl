using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("gm_complexallostery.jl"))

################################################################################
# The simple case, C=2 and B mostly 1
################################################################################
function make_v2(N, B; edge_t=:full)
    ca = ComplexAllosteryGM(N, 2, B;
        energy_matrices=make_EM_sym_C2(B),
        version=2
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
get_rs() = Symbolics.variables(:r, 1:3)
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

                E_state = calc_energy(state, ca)
                E_new_state = calc_energy(new_state, ca)

                # forward rate
                exp_term_f = thetas[3, 1] * E_state - (1 - thetas[3, 2]) * E_new_state
                rf = rs[3] * exp(exp_term_f / kT)
                add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                # backward rate
                exp_term_b = thetas[3, 2] * E_new_state - (1 - thetas[3, 1]) * E_state
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

                E_state = calc_energy(state, ca)
                E_new_state = calc_energy(new_state, ca)

                # P + S <-> NS reaction
                begin
                    # forward rate
                    exp_term_f = thetas[1, 1] * (E_state + chem_energies[1]) - (1 - thetas[1, 2]) * E_new_state
                    rf = rs[1] * concetrations[1] * exp(exp_term_f / kT)
                    add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[1, 2] * E_new_state - (1 - thetas[1, 1]) * (E_state + chem_energies[1])
                    rb = rs[1] * exp(exp_term_b / kT)
                    add_edge_incweight!(ca.graph, new_vertex, vertex, rb)
                end

                # ATP + S <-> ADP + NS reaction
                begin
                    # forward rate
                    exp_term_f = thetas[2, 1] * (E_state + chem_energies[2]) - (1 - thetas[2, 2]) * (E_new_state + chem_energies[3])
                    rf = rs[2] * concetrations[2] * exp(exp_term_f / kT)
                    add_edge_incweight!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[2, 2] * (E_new_state + chem_energies[3]) - (1 - thetas[2, 1]) * (E_state + chem_energies[1])
                    rb = rs[2] * concetrations[3] * exp(exp_term_b / kT)
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
    terms

    substitute_partial(ca, terms)
end
