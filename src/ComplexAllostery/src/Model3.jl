module Model3
using ..ComplexAllostery

################################################################################
# Setting up the variables
################################################################################
# These variables define the system
struct SystemVars{F<:Number}
    energies::SVector{3,F} # ε_t, Δε_r, ε_b in order
    rs::Tuple{SVector{2,F},SVector{2,F},F}
    chem_energies::SVector{3,F} # in order of P, ATP, ADP
    concentrations::SVector{3,F}
    thetas::SMatrix{3,2,F} # first means forwards, second backwards
    kT::F
end
function smart_vars(
    energies::SVector{3,A},
    rs::Tuple{SVector{2,B},SVector{2,C},D},
    chem_energies::SVector{3,F},
    concentrations::SVector{3,E},
    thetas::SMatrix{3,2,G},
    kT::H;
    concretetype::Val{FT}=Val(Float64)
) where {A,B,C,D,E,F,G,H,FT}
    if Num in SA[A, B, C, D, E, F, G, H]
        SystemVars(
            convert.(Num, energies),
            (convert.(Num, rs[1]), convert.(Num, rs[2]), convert(Num, rs[3])),
            convert.(Num, chem_energies),
            convert.(Num, concentrations),
            convert.(Num, thetas),
            convert(Num, kT)
        )
    else
        SystemVars(
            convert.(FT, energies),
            (convert.(FT, rs[1]), convert.(FT, rs[2]), convert(FT, rs[3])),
            convert.(FT, chem_energies),
            convert.(FT, concentrations),
            convert.(FT, thetas),
            convert(FT, kT)
        )
    end
end
export SystemVars, smart_vars

# The major symbolic variables so that I don't have to keep redefining them
symvars = (;
    # The base model 2.5 variables
    energies=listsavariables(:ε_t, :Δε_r, :ε_b),
    rs=(savariables(:r, 1, 1:2), savariables(:r, 2, 1:2), variable(:r, 3)),
    chem_energies=listsavariables(:εP, :εATP, :εADP),
    concentrations=listsavariables(:cP, :cATP, :cADP),
    thetas=savariables(:θ, 1:3, 1:2),
    # Simplifed model vars
    sim_rs=SVector(variables(:r, 1:3)..., variable(:α))
)
export symvars

function vars_base()
    smart_vars(
        symvars.energies,
        symvars.rs,
        symvars.chem_energies,
        symvars.concentrations,
        symvars.thetas,
        symvar_kT
    )
end
export vars_base

function vars_simplified(;
    newrs=:basic,
    energies=:onlyeb,
    concentrations=nothing,
    kT=true
)
    base = vars_base()

    if isnothing(newrs)
        newrs = base.rs
    else
        if newrs == :basic
            rs = symvars.sim_rs[1:3]
            alpha = symvars.sim_rs[4]
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
    end

    if isnothing(energies)
        energies = base.energies
    elseif energies == :onlyeb
        energies = SA[0.0, 0.0, base.energies[end]]
    elseif ((isa(energies, AbstractVector) || isa(energies, Tuple))) && (length(energies) == 2)
        energies = SVector(energies[1], energies[2], base.energies[end])
    else
        throw(ArgumentError(f"invalid energies of {energies}"))
    end

    if isnothing(concentrations)
        concentrations = base.concentrations
    elseif concentrations == :ss1 # cADP -> 0
        concentrations = setindex(base.concentrations, 0.0, 3)
    elseif concentrations == :ss2 # cADP -> 0, cP -> 0
        concentrations = setindex(base.concentrations, 0.0, 3)
        concentrations = setindex(concentrations, 0.0, 1)
    else
        throw(ArgumentError(f"invalid concentrations of {concentrations}"))
    end

    if kT == true
        kT = 1.0
    elseif kT == false
        kT = symvar_kT
    else
        throw(ArgumentError(f"invalid kT of {kT}"))
    end

    smart_vars(
        energies,
        newrs,
        (@SVector fill(Num(0.0), 3)),
        concentrations,
        (@SMatrix fill(Num(1.0), 3, 2)),
        kT
    )
end
export vars_simplified

################################################################################
# Making the system given SystemVars
################################################################################
function makeEM(et::F, der::F, eb::F, B=1) where {F}
    monomer = Matrix{F}(undef, 2, B + 1)
    monomer[1, 1] = 0
    monomer[2, 1] = der
    monomer[1, 2:B+1] .= et .* collect(1:B)
    monomer[2, 2:B+1] .= monomer[1, 2:B+1] .- der

    interactions = [0 eb; eb 0]
    EnergyMatrices(monomer, interactions)
end
function makeEM(vars::SystemVars{F}=vars_base(), B=1) where {F}
    makeEM(vars.energies[1], vars.energies[2], vars.energies[3], B)
end

function makeCAGM(N, vars::SystemVars{F}=vars_base(), symmetry=Loop(), B=1; kwargs...) where {F}
    ca = ComplexAllosteryGM(N, 2, B;
        symmetry,
        numtype=Val(F),
        energy_matrices=makeEM(vars, B),
        version=3.0,
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

                almost_rf = exp((vars.thetas[3, 1] * e_i_state - (1 - vars.thetas[3, 2]) * e_i_new_state) / vars.kT)
                almost_rb = exp((vars.thetas[3, 2] * e_i_new_state - (1 - vars.thetas[3, 1]) * e_i_state) / vars.kT)

                divisor = almost_rf + almost_rb

                rf = vars.rs[3] * almost_rf / divisor
                inc_edge!(ca.graph, vertex, new_vertex, rf)
                rb = vars.rs[3] * almost_rb / divisor
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
                    almost_rf = vars.concentrations[1] * exp(exp_term_f / vars.kT)
                    # backward rate
                    exp_term_b = vars.thetas[1, 2] * e_i_new_state - (1 - vars.thetas[1, 1]) * (e_i_state + vars.chem_energies[1])
                    almost_rb = exp(exp_term_b / vars.kT)

                    divisor = almost_rf + almost_rb

                    rf = vars.rs[1][i_conf] * almost_rf / divisor
                    inc_edge!(ca.graph, vertex, new_vertex, rf)
                    rb = vars.rs[1][i_conf] * almost_rb / divisor
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end

                # ATP + S <-> ADP + NS reaction
                begin
                    # forward rate
                    exp_term_f = vars.thetas[2, 1] * (e_i_state + vars.chem_energies[2]) - (1 - vars.thetas[2, 2]) * (e_i_new_state + vars.chem_energies[3])
                    almost_rf = vars.concentrations[2] * exp(exp_term_f / vars.kT)
                    # backward rate
                    exp_term_b = vars.thetas[2, 2] * (e_i_new_state + vars.chem_energies[3]) - (1 - vars.thetas[2, 1]) * (e_i_state + vars.chem_energies[2])
                    almost_rb = vars.concentrations[3] * exp(exp_term_b / vars.kT)

                    divisor = almost_rf + almost_rb

                    rf = vars.rs[2][i_conf] * almost_rf / divisor
                    inc_edge!(ca.graph, vertex, new_vertex, rf)
                    rb = vars.rs[2][i_conf] * almost_rb / divisor
                    inc_edge!(ca.graph, new_vertex, vertex, rb)
                end
            end
        end
    end

    ca
end
export makeCAGM

################################################################################
# Helper functions for later substitutions
################################################################################
function terms_simplified(basevars=vars_base(); kwargs...)
    simvars = vars_simplified(; kwargs...)

    kaka1 = [basevars.energies; basevars.rs[1]; basevars.rs[2]; basevars.rs[3]; basevars.chem_energies; basevars.concentrations; collect(Iterators.flatten(basevars.thetas)); basevars.kT]
    kaka2 = [simvars.energies; simvars.rs[1]; simvars.rs[2]; simvars.rs[3]; simvars.chem_energies; simvars.concentrations; collect(Iterators.flatten(simvars.thetas)); simvars.kT]

    terms = Dict{Num,Num}()
    for (old, new) in zip(kaka1, kaka2)
        if !isequal(old, new)
            terms[old] = new
        end
    end
    terms
end
export terms_simplified

function terms_delta_concentrations()
    conc = vars_base().concentrations
    dconc = Symbolics.variable.([:dP, :dATP, :dADP])

    terms = Dict{Num,Num}()
    for (c, dc) in zip(conc, dconc)
        terms[c] = 1 + dc
    end

    terms
end
export terms_delta_concentrations

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
export terms_equilibrium

function terms_energies(; et=nothing, der=nothing, eb=nothing)
    basevars = vars_base()
    terms = Dict{Num,Num}()
    if !isnothing(et)
        terms[basevars.energies[1]] = et
    end
    if !isnothing(der)
        terms[basevars.energies[2]] = der
    end
    if !isnothing(eb)
        terms[basevars.energies[3]] = eb
    end
    terms
end
terms_energies(et, der=nothing, eb=nothing) = terms_energies(; et, der, eb)
export terms_energies

end
