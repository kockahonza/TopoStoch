using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gmca.jl"))

################################################################################
# The C=2 case for Model 2.5
################################################################################
struct Vars_2_5_C2{F<:Number}
    energies::SVector{3,F} # ε_t, Δε_r, ε_b in order
    rs::Tuple{SVector{2,F},SVector{2,F},F}
    chem_energies::SVector{3,F} # in order of P, ATP, ADP
    concentrations::SVector{3,F}
    thetas::SMatrix{3,2,F} # first means forwards, second backwards
    kT::F
end
function make_vars_smart(
    energies::SVector{3,A},
    rs::Tuple{SVector{2,B},SVector{2,C},D},
    chem_energies::SVector{3,F},
    concentrations::SVector{3,E},
    thetas::SMatrix{3,2,G},
    kT::H;
    concretetype::Val{FT}=Val(Float64)
) where {A,B,C,D,E,F,G,H,FT}
    if Num in SA[A, B, C, D, E, F, G, H]
        Vars_2_5_C2(
            convert.(Num, energies),
            (convert.(Num, rs[1]), convert.(Num, rs[2]), convert(Num, rs[3])),
            convert.(Num, chem_energies),
            convert.(Num, concentrations),
            convert.(Num, thetas),
            convert(Num, kT)
        )
    else
        Vars_2_5_C2(
            convert.(FT, energies),
            (convert.(FT, rs[1]), convert.(FT, rs[2]), convert(FT, rs[3])),
            convert.(FT, chem_energies),
            convert.(FT, concentrations),
            convert.(FT, thetas),
            convert(FT, kT)
        )
    end
end
function make_sem_C2(B, vars::Vars_2_5_C2)
    make_sem_C2(B; et=vars.energies[1], der=vars.energies[2], eb=vars.energies[3])
end

################################################################################
# Setting up particular Vars_2_5_C2 choices
################################################################################
# Making simplified models directly using Vars_2_5_C2
function vars_base()
    r12s = Symbolics.variables(:r, 1:2, 1:2)
    rs = [SVector(r12s[i, :]...) for i in 1:2]
    make_vars_smart(
        get_sem_C2_vars(),
        (rs[1], rs[2], Symbolics.variable(:r, 3)),
        SVector(Symbolics.variable.([:εP, :εATP, :εADP])...),
        SVector(Symbolics.variable.([:cP, :cATP, :cADP])...),
        SMatrix{3,2}(Symbolics.variables(:θ, 1:3, 1:2)...),
        Symbolics.variable(:kT)
    )
end

function vars_simplified(;
    newrs=nothing,
    energies=nothing,
    concentrations=nothing,
    kT=nothing
)
    if isnothing(newrs)
        rs = Symbolics.variables(:r, 1:3)
        alpha = Symbolics.variable(:α)
    elseif newrs == :ss1
        rs = [1.0, 1.0, 1.0]
        alpha = 0.0
    elseif ((isa(newrs, AbstractVector) || isa(newrs, Tuple)) && (length(newrs) == 4))
        rs = newrs[1:3]
        alpha = newrs[4]
    else
        throw(ArgumentError(f"invalid newrs of {newrs}"))
    end
    newrs = (SVector(rs[1], rs[1]), SVector(rs[2], alpha * rs[2]), rs[3])

    if isnothing(energies) || (energies == false)
        energies = get_sem_C2_vars()
    elseif energies == :noall
        energies = SA[0.0, 0.0, get_sem_C2_vars()[end]]
    elseif (isa(energies, AbstractVector) || isa(energies, Tuple))
        if length(energies) == 2
            energies = SVector(energies[1], energies[2], get_sem_C2_vars()[end])
        end
    else
        throw(ArgumentError(f"invalid energies of {energies}"))
    end

    if isnothing(concentrations) || (concentrations == false)
        concentrations = SVector(Symbolics.variable.([:cP, :cATP, :cADP])...)
    elseif concentrations == :ss1
        concentrations = SVector(Symbolics.variable(:cP), Symbolics.variable(:cATP), 0.0)
    elseif concentrations == :ss2
        concentrations = SVector(0.0, Symbolics.variable(:cATP), 0.0)
    else
        throw(ArgumentError(f"invalid concentrations of {concentrations}"))
    end

    if isnothing(kT)
        kT = 1.0
    elseif (kT == :sym) || (kT == false)
        kT = Symbolics.variable(:kT)
    else
        throw(ArgumentError(f"invalid kT of {kT}"))
    end

    make_vars_smart(
        energies,
        newrs,
        (@SVector fill(Num(0.0), 3)),
        concentrations,
        (@SMatrix fill(Num(1.0), 3, 2)),
        kT
    )
end

################################################################################
# Making the ComplexAllosteryGM
################################################################################
function make_v2_5(N, B;
    symmetry=Loop(),
    vars::Vars_2_5_C2{F}=vars_base(),
    kwargs...
) where {F}
    ca = ComplexAllosteryGM(N, 2, B;
        symmetry,
        energy_matrices=make_sem_C2(B, vars),
        version=2.5,
        metadata=Dict("vars" => vars),
        kwargs...
    )

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
                exp_term_f = vars.thetas[3, 1] * e_i_state - (1 - vars.thetas[3, 2]) * e_i_new_state
                rf = vars.rs[3] * exp(exp_term_f / vars.kT)
                inc_edge!(ca.graph, vertex, new_vertex, rf)

                # backward rate (where c_i decreases)
                exp_term_b = vars.thetas[3, 2] * e_i_new_state - (1 - vars.thetas[3, 1]) * e_i_state
                rb = vars.rs[3] * exp(exp_term_b / vars.kT)
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
                    exp_term_f = vars.thetas[1, 1] * (e_i_state + vars.chem_energies[1]) - (1 - vars.thetas[1, 2]) * e_i_new_state
                    rf = vars.rs[1][i_conf] * vars.concentrations[1] * exp(exp_term_f / vars.kT)
                    inc_edge!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = vars.thetas[1, 2] * e_i_new_state - (1 - vars.thetas[1, 1]) * (e_i_state + vars.chem_energies[1])
                    rb = vars.rs[1][i_conf] * exp(exp_term_b / vars.kT)
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end

                # ATP + S <-> ADP + NS reaction
                begin
                    # forward rate
                    exp_term_f = vars.thetas[2, 1] * (e_i_state + vars.chem_energies[2]) - (1 - vars.thetas[2, 2]) * (e_i_new_state + vars.chem_energies[3])
                    rf = vars.rs[2][i_conf] * vars.concentrations[2] * exp(exp_term_f / vars.kT)
                    inc_edge!(ca.graph, vertex, new_vertex, rf)

                    # backward rate
                    exp_term_b = vars.thetas[2, 2] * (e_i_new_state + vars.chem_energies[3]) - (1 - vars.thetas[2, 1]) * (e_i_state + vars.chem_energies[2])
                    rb = vars.rs[2][i_conf] * vars.concentrations[3] * exp(exp_term_b / vars.kT)
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end
            end
        end
    end

    ca
end

# Also make terms for further simplifying using ssubstitute and similar
function terms_simplified(; kwargs...)
    basevars = vars_base()
    simvars = vars_simplified(; kwargs...)

    oldterms = Num[]
    newterms = Num[]

    if !isequal(basevars.energies, simvars.energies)
        append!(oldterms, basevars.energies)
        append!(newterms, simvars.energies)
    end

    for i in 1:3
        append!(oldterms, basevars.rs[i])
        append!(newterms, simvars.rs[i])
    end
    append!(oldterms, basevars.chem_energies)
    append!(newterms, simvars.chem_energies)
    if !isequal(basevars.concentrations, simvars.concentrations)
        append!(oldterms, basevars.concentrations)
        append!(newterms, simvars.concentrations)
    end
    append!(oldterms, basevars.thetas)
    append!(newterms, simvars.thetas)
    if !isequal(basevars.kT, simvars.kT)
        append!(oldterms, basevars.kT)
        append!(newterms, simvars.kT)
    end

    terms = Dict{Num,Num}()
    for (old, new) in zip(oldterms, newterms)
        terms[old] = new
    end
    terms
end

function terms_delta_concentrations()
    conc = vars_base().concentrations
    dconc = Symbolics.variable.([:dP, :dATP, :dADP])

    terms = Dict{Num,Num}()
    for (c, dc) in zip(conc, dconc)
        terms[c] = 1 + dc
    end

    terms
end

function terms_equilibrium(elim=:cATP)
    terms = Dict{Num,Num}()
    basevars = vars_base()
    es = basevars.chem_energies
    cs = basevars.concentrations
    if elim == :cATP
        terms[es[2]] = es[1] + es[3]
        terms[cs[2]] = cs[1] * cs[3]
    elseif elim == :cP
        terms[es[1]] = es[2] - es[3]
        terms[cs[1]] = cs[2] / cs[3]
    elseif elim == :cADP
        terms[es[2]] = es[1] - es[3]
        terms[cs[2]] = cs[1] / cs[3]
    end
    terms
end
