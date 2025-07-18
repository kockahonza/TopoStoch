@reexport using MetaGraphsNext

struct MolAut{S<:Symmetry} <: AbstractGraphModel{Float64}
    L::Int
    reduced::Bool
    mg::MetaGraph{Int,SimpleDiGraph{Int},Vector{Int},String,Float64}
end
function MolAut(L, rule::Int; reduced=false)
    ned = make_ca_ned(L, rule; show=false)
    g = graph(ned)
    states = collect(allstates(ned))

    mg = MetaGraph(
        DiGraph();
        label_type=Vector{Int},
        vertex_data_type=String,
        edge_data_type=Float64,
        weight_function=identity,
        default_weight=0.0,
        graph_data=(;
            rule, L,
            desc="all state ned graph for rule $rule and L=$L"
        )
    )

    for s in states
        mg[s] = join(s)
    end

    for e in edges(g)
        mg[states[e.src], states[e.dst]] = get_weight(g, e.src, e.dst)
    end

    if reduced
        mg = group_shifted_states_mg(mg)
    end

    mg

    MolAut{Loop}(L, reduced, mg)
end
function graph(ma::MolAut)
    swg = SimpleWeightedDiGraph{Int,Float64}(nv(ma.mg))
    for e in edges(ma.mg)
        add_edge!(swg, e.src, e.dst, ma.mg[label_for(ma.mg, e.src), label_for(ma.mg, e.dst)])
    end
    swg
end
numstates(ma::MolAut) = nv(ma.mg)
allstates(ma::MolAut) = collect(labels(ma.mg))
export MolAut

################################################################################
# Reducing a ned/ma graph
################################################################################
function group_shifted_states_mg(mg)
    nmg_metadata = mg[]

    grouped_nmg = MetaGraph(
        DiGraph();
        label_type=Vector{Int},
        vertex_data_type=String,
        edge_data_type=Float64,
        weight_function=identity,
        default_weight=0.0,
        graph_data=(;
            mg[]..., desc="merged ned graph for rule $(nmg_metadata.rule) and L=$(nmg_metadata.L)"
        )
    )

    states = labels(mg)
    ustates, state_to_ustate_is = arg_group_ned_states(states)

    for us in ustates
        grouped_nmg[us] = join(us)
    end

    for e in edges(mg)
        sc = e.src
        dc = e.dst
        s = label_for(mg, sc)
        d = label_for(mg, dc)
        us = ustates[state_to_ustate_is[sc]]
        ud = ustates[state_to_ustate_is[dc]]
        if !haskey(grouped_nmg, us, ud)
            grouped_nmg[us, ud] = mg[s, d]
        else
            grouped_nmg[us, ud] += mg[s, d]
        end
    end

    grouped_nmg
end

function states_related_by_shift(st1, st2)
    if length(st1) != length(st2)
        return false
    elseif sum(st1) != sum(st2)
        return false
    else
        for o in 0:length(st1)-1
            are_equal = true
            for i in 1:length(st1)
                if st1[i] != st2[mod1(i + o, length(st2))]
                    are_equal = false
                    break
                end
            end
            if are_equal
                return true
            end
        end
        return false
    end
end

function arg_group_ned_states(states)
    unique_rules = []
    group_is = fill(-1, length(states))
    for (st_i, st) in enumerate(states)
        found_previous = false
        for (ur_i, ur) in enumerate(unique_rules)
            if states_related_by_shift(st, ur)
                found_previous = true
                group_is[st_i] = ur_i
                if st > ur
                    unique_rules[ur_i] = st
                end
                break
            end
        end
        if !found_previous
            push!(unique_rules, st)
            group_is[st_i] = length(unique_rules)
        end
    end
    unique_rules, group_is
end
