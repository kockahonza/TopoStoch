using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gmca.jl"))

################################################################################
# The C=2 case for Model 2.5
################################################################################
# define a bunch of symbolic variables so I don't have to keep redoing this
function get_rs()
    r12s = Symbolics.variables(:r, 1:2, 1:2)
    rs = [r12s[i, :] for i in 1:2]
    [rs; Symbolics.variable(:r, 3)]
end
get_thetas() = Symbolics.variables(:θ, 1:3, 1:2)
get_kT() = Symbolics.variable(:kT)
get_chem_energies() = Symbolics.variable.([:εP, :εATP, :εADP])
get_concentrations() = Symbolics.variable.([:cP, :cATP, :cADP])

function make_v2_5(N, B; symmetry=Loop(), simplified=false, noall=false, kwargs...)
    ca = ComplexAllosteryGM(N, 2, B;
        symmetry,
        energy_matrices=noall ? make_sem_C2_noall(B) : make_sem_C2(B),
        version=2.5, kwargs...
    )

    # Do the edges
    # Declare symbolic variables here
    rs = get_rs()
    thetas = get_thetas()
    kT = get_kT()
    chem_energies = get_chem_energies()
    concentration = get_concentrations()

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
                inc_edge!(ca.graph, vertex, new_vertex, rf)

                # backward rate (where c_i decreases)
                exp_term_b = thetas[3, 2] * e_i_new_state - (1 - thetas[3, 1]) * e_i_state
                rb = rs[3] * exp(exp_term_b / kT)
                inc_edge!(ca.graph, new_vertex, vertex, rb)
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
                    rf = rs[1][i_conf] * concentration[1] * exp(exp_term_f / kT)
                    inc_edge!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[1, 2] * e_i_new_state - (1 - thetas[1, 1]) * (e_i_state + chem_energies[1])
                    rb = rs[1][i_conf] * exp(exp_term_b / kT)
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end

                # ATP + S <-> ADP + NS reaction
                begin
                    # forward rate
                    exp_term_f = thetas[2, 1] * (e_i_state + chem_energies[2]) - (1 - thetas[2, 2]) * (e_i_new_state + chem_energies[3])
                    rf = rs[2][i_conf] * concentration[2] * exp(exp_term_f / kT)
                    inc_edge!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = thetas[2, 2] * (e_i_new_state + chem_energies[3]) - (1 - thetas[2, 1]) * (e_i_state + chem_energies[2])
                    rb = rs[2][i_conf] * concentration[3] * exp(exp_term_b / kT)
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end
            end
        end
    end

    if simplified == :both
        ca, r_subs(ca)
    elseif simplified == :full
        r_subs(ca; r1=1.0, r2=1.0, r3=1.0, alpha=0.0)
    elseif simplified
        r_subs(ca)
    elseif !simplified
        ca
    end
end

################################################################################
# Simplifications and running
################################################################################
"""
Simplifies by setting all thetas to 1 and chemical energies to 0.
"""
function chem_subs(obj; kT=true)
    terms = Dict{Num,Num}()
    for theta in get_thetas()
        terms[theta] = 1.0
    end
    for chem_energy in get_chem_energies()
        terms[chem_energy] = 0.0
    end

    if kT
        terms[get_kT()] = 1.0
    end

    ssubstitute(obj, terms)
end

"""
Simplified by setting both r1s to r1, r2(1) = r2 and r2(2)=r2*alpha,
can also do chem_subs.
"""
function get_newrs()
    nrs = Symbolics.variables(:r, 1:3)
    [nrs; Symbolics.variable(:α)]
end
function r_subs(obj; do_chem_subs=true, r1=nothing, r2=nothing, r3=nothing, alpha=nothing, kwargs...)
    rs = get_rs()
    newrs = Symbolics.variables(:r, 1:3)
    if isnothing(r1)
        r1 = newrs[1]
    end
    if isnothing(r2)
        r2 = newrs[2]
    end
    if isnothing(r3)
        r3 = newrs[3]
    end
    if isnothing(alpha)
        alpha = Symbolics.variable(:α)
    end

    terms = Dict{Num,Num}()

    terms[rs[1][1]] = r1
    terms[rs[1][2]] = r1

    terms[rs[2][1]] = r2
    terms[rs[2][2]] = alpha * r2

    terms[rs[3]] = r3

    if do_chem_subs
        obj = chem_subs(obj; kwargs...)
    end
    ssubstitute(obj, terms)
end
function ssr_terms()
    terms = Dict(get_newrs()[1:3] .=> 1.0)
    terms[get_newrs()[4]] = 0.0
    terms
end

function dconc_terms()
    conc = get_concentrations()
    dconc = Symbolics.variable.([:dP, :dATP, :dADP])

    terms = Dict{Num,Num}()
    for (c, dc) in zip(conc, dconc)
        terms[c] = 1 + dc
    end

    terms
end

function int_plot_1(ca::ComplexAllosteryGM, args...; init_scen=:eonly, kwargs...)
    if ca.version != 2.5
        throw(ArgumentError("this method only works for version 2.5"))
    end
    sca, frterms = r_subs(ca)
    variables = vcat(get_em_vars(), get_concentrations(), frterms)
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
    terms[es[2]] = es[1] + es[3]
    cs = get_concentrations()
    terms[cs[2]] = cs[1] * cs[3]
    terms
end
