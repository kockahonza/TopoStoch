################################################################################
# Dealing with energy and equilibrium probability calculations
################################################################################
# Two symbolics terms that might be needed throughout
const symvar_mu = Symbolics.variable(:μ)
const symvar_kT = Symbolics.variable(:kT)
export symvar_mu, symvar_kT

# Simple calculator functions
function get_neighbors(st::CAState, i, ::Chain)
    N = length(st.conformations)
    if N == 1
        return (nothing, nothing)
    end

    if i == 1
        (st.conformations[i+1], nothing)
    elseif i == N
        (st.conformations[N-1], nothing)
    else
        (st.conformations[i-1], st.conformations[i+1])
    end
end
function get_neighbors(st::CAState, i, ::Loop)
    N = length(st.conformations)
    if N == 1
        return (nothing, nothing)
    end

    if i == 1
        (st.conformations[N], st.conformations[i+1],)
    elseif i == N
        (st.conformations[N-1], st.conformations[1])
    else
        (st.conformations[i-1], st.conformations[i+1])
    end
end
function get_neighbors(st::CAState, i, _::ComplexAllosteryGM{S}) where {S}
    get_neighbors(st, i, S())
end
export get_neighbors

function calc_monomer_energy(st::CAState, i, em::EnergyMatrices)
    em.monomer[st.conformations[i], st.occupations[i]+1]
end
calc_monomer_energy(st::CAState, i, ca::ComplexAllosteryGM) = calc_monomer_energy(st, i, ca.energy_matrices)
export calc_monomer_energy

function calc_interaction_energy(st::CAState, i, em::EnergyMatrices, s::Symmetry)
    (left, right) = get_neighbors(st, i, s)
    energy = 0
    if !isnothing(left)
        energy += em.interactions[left, st.conformations[i]]
    end
    if !isnothing(right)
        energy += em.interactions[st.conformations[i], right]
    end
    energy
end
calc_interaction_energy(st::CAState, i, ca::ComplexAllosteryGM{S}) where {S} = calc_interaction_energy(st, i, ca.energy_matrices, S())
export calc_interaction_energy

function calc_energy(st::CAState, em::EnergyMatrices, s::Symmetry)
    energy = 0
    for i in 1:length(st.conformations)
        energy += calc_monomer_energy(st, i, em)
        energy += 0.5 * calc_interaction_energy(st, i, em, s)
    end
    energy
end
function calc_energy(st::CAState, ca::ComplexAllosteryGM{S}) where {S}
    calc_energy(st, ca.energy_matrices, S())
end
export calc_energy

calc_numligands(st::CAState) = sum(st.occupations)
export calc_numligands

calc_numofconf(st::CAState, conf_state=1) = count(x -> x == conf_state, st.conformations)
export calc_numofconf

function calc_numboundaries(st::CAState, args...)
    total = 0.0
    for i in 1:length(st.conformations)
        neighbors = filter(!isnothing, get_neighbors(st, i, args...))
        total += 0.5 * count(x -> x != st.conformations[i], neighbors)
    end
    Int(total)
end
export calc_numboundaries

function calc_gibbs_factor(st::CAState, em::EnergyMatrices,
    s::Symmetry, μ=symvar_mu, kT=symvar_kT)
    exp(-(calc_energy(st, em, s) - μ * calc_numligands(st)) / kT)
end
function calc_gibbs_factor(st::CAState, ca::ComplexAllosteryGM{S}, args...) where {S}
    env = if length(args) > 0
        args
    elseif !isnothing(ca.environment)
        ca.environment
    else
        ()
    end
    calc_gibbs_factor(st, ca.energy_matrices, S(), env...)
end
export calc_gibbs_factor

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
export calc_all_gibbs_factors, calc_partition_function, calc_avg_numligands, calc_avg_energy

################################################################################
# Dealing with transitions
################################################################################
# Adding some basic links without meaningful weights
function add_edges_base!(ca::ComplexAllosteryGM)
    # Just use a constant for all the rates, disregarding DB and all else too
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

function add_edges_balanced_universtal_scale!(ca::ComplexAllosteryGM)
    universal_scale = Symbolics.variable(:r)
    gibbs_factors = calc_gibbs_factor.(allstates(ca), ca)

    function add_new_rate(v1, v2)
        if (v1 < v2)
            divp = gibbs_factors[v2] / gibbs_factors[v1]
            expfactor = sqrt(divp)
            add_edge!(ca.graph, v1, v2, simplify(universal_scale * expfactor))
            add_edge!(ca.graph, v2, v1, simplify(universal_scale / expfactor))
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

function add_edges_fcs_energetics!(ca::ComplexAllosteryGM; weight=1)
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

                if e_i_new_state < e_i_state
                    inc_edge!(ca.graph, vertex, new_vertex, weight)
                else
                    inc_edge!(ca.graph, new_vertex, vertex, weight)
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                new_vertex = statetoi(new_state, ca)

                i_conf = state.conformations[i]

                if i_conf == 1
                    inc_edge!(ca.graph, vertex, new_vertex, weight)
                else
                    inc_edge!(ca.graph, new_vertex, vertex, weight)
                end
            end
        end
    end
end

function add_edges_fcs_simple!(ca::ComplexAllosteryGM; weight=1)
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)

        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in (state.conformations[i]+1):ca.C
                new_state = copy(state)
                new_state.conformations[i] = new_c
                new_vertex = statetoi(new_state, ca)

                bi = state.occupations[i]

                if bi == 0
                    inc_edge!(ca.graph, new_vertex, vertex, weight)
                else
                    inc_edge!(ca.graph, vertex, new_vertex, weight)
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                new_vertex = statetoi(new_state, ca)

                i_conf = state.conformations[i]

                if i_conf == 1
                    inc_edge!(ca.graph, vertex, new_vertex, weight)
                else
                    inc_edge!(ca.graph, new_vertex, vertex, weight)
                end
            end
        end
    end
end

function make_simple_fcs_only(N, B; edge_t=:simple, kwargs...)
    ca = ComplexAllosteryGM(N, 2, B;
        symmetry=Loop(),
        energy_matrices=nothing,
        version=42.0,
        numtype=Val(Float64),
        kwargs...
    )

    if edge_t == :simple
        add_edges_fcs_simple!(ca)
    elseif edge_t == :en
        add_edges_fcs_energetics!(ca)
    else
        throw(ArgumentError(f"edge_t \"{edge_t}\" not recognised"))
    end

    ca
end

"""
Check detailed balance using equilibrium stats which will include μ
"""
function check_eq_db(ca::ComplexAllosteryGM; sim=false, fil=true)
    gfs = calc_gibbs_factor.(allstates(ca), ca)

    conditions = []
    for i in 1:numstates(ca)
        for j in i:numstates(ca)
            wij = get_weight(ca.graph, i, j)
            wji = get_weight(ca.graph, j, i)
            zero1 = iszero(wij)
            zero2 = iszero(wji)
            if !zero1 && !zero2
                push!(conditions, (gfs[i] * wij) - (gfs[j] * wji))
            elseif !(zero1 && zero2)
                @error f"detailed balance is not satisfied as some non-zero transitions do not have a reverse transition! (i, j)={(i, j)}"
            end
        end
    end

    if sim
        for i in eachindex(conditions)
            kaka = simplify(conditions[i]; rewriter=get_rewriter())
            if !iszero(kaka)
                println(kaka)
            end
        end
    end

    if fil
        conditions = filter(!iszero, conditions)
    end

    conditions
end
