using DrWatson;
@quickactivate "TopoStochSim";

using GraphModels
using NonEqDigits
using ColorSchemes
using GraphvizDotLang: digraph, subgraph, node, edge
using MetaGraphsNext

function find_singlestateac_rules(L, codes=ca_ucodes_f1())
    filter(codes) do c
        ned = make_ca_ned(L, c; show=false)
        g = graph(ned)
        all(attracting_components(g)) do ac
            length(ac) == 1
        end
    end
end

function normalized_adjacency_matrix(g; mat=true)
    M = adjacency_matrix(g)
    if mat
        M = collect(M)
    end
    for i in 1:nv(g)
        Mrow = @view M[i, :]
        ss = sum(Mrow)
        if !iszero(ss)
            Mrow ./= ss
        end
    end
    M
end
normalized_adjacency_matrix(gm::AbstractGraphModel) = normalized_adjacency_matrix(graph(gm))

################################################################################
# Version 2.5 - accounting for shifts/rotations
################################################################################
function make_ned_meta_graph(L, code;
    merge_shifted=false
)
    ned = make_ca_ned(L, code; show=false)
    g = graph(ned)
    states = collect(allstates(ned))

    nmg = MetaGraph(
        DiGraph();
        label_type=Vector{Int},
        vertex_data_type=String,
        edge_data_type=Float64,
        weight_function=identity,
        default_weight=0.0,
        graph_data=(;
            code, L,
            desc="all state ned graph for code $code and L=$L"
        )
    )

    for s in states
        nmg[s] = join(s)
    end

    for e in edges(g)
        nmg[states[e.src], states[e.dst]] = get_weight(g, e.src, e.dst)
    end

    if merge_shifted
        group_shifted_states_mg(nmg)
    else
        nmg
    end
end

function group_shifted_states_mg(nmg)
    nmg_metadata = nmg[]

    grouped_nmg = MetaGraph(
        DiGraph();
        label_type=Vector{Int},
        vertex_data_type=String,
        edge_data_type=Float64,
        weight_function=identity,
        default_weight=0.0,
        graph_data=(;
            nmg[]..., desc="merged ned graph for code $(nmg_metadata.code) and L=$(nmg_metadata.L)"
        )
    )

    states = labels(nmg)
    ustates, state_to_ustate_is = arg_group_ned_states(states)

    for us in ustates
        grouped_nmg[us] = join(us)
    end

    for e in edges(nmg)
        sc = e.src
        dc = e.dst
        s = label_for(nmg, sc)
        d = label_for(nmg, dc)
        us = ustates[state_to_ustate_is[sc]]
        ud = ustates[state_to_ustate_is[dc]]
        if !haskey(grouped_nmg, us, ud)
            grouped_nmg[us, ud] = nmg[s, d]
        else
            grouped_nmg[us, ud] += nmg[s, d]
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

function group_ned_states(states)
    results = Dict{eltype(states),Vector{eltype(states)}}()
    for st in states
        found_previous = false
        for k in keys(results)
            if states_related_by_shift(st, k)
                found_previous = true
                push!(results[k], st)
                break
            end
        end
        if !found_previous
            results[st] = [st]
        end
    end
    collect(values(results))
end

function plot_nmg!(ax, nmg;
    nlabels=true,
    ecolors=true,
    edge_attr=(;),
    kwargs...
)
    auto_kwargs = Dict{Symbol,Any}()
    if nlabels
        auto_kwargs[:nlabels] = [nmg[l] for l in labels(nmg)]
    end
    if ecolors && !haskey(kwargs, :edge_color)
        edge_color = [nmg[s, d] for (s, d) in edge_labels(nmg)]
        ecmin, ecmax = extrema(edge_color)
        if !(ecmin == ecmax)
            auto_kwargs[:edge_color] = edge_color
        else
            @warn "All edges have the same value, skipping colors"
        end
    end

    auto_kwargs[:edge_attr] = edge_attr
    auto_kwargs[:arrow_attr] = edge_attr

    graphplot!(ax, nmg; auto_kwargs..., kwargs...)
end

function plot_nmg(nmg;
    layout=Spring(dim=2),
    add_colorbars=true,
    kwargs...
)
    if !isa(layout, AbstractArray)
        layout = layout(nmg)
    end
    dim = length(layout[1])

    fig = Figure()
    ax = if dim == 2
        Axis(fig[1, 1])
    else
        LScene(fig[1, 1])
    end
    p = plot_nmg!(ax, nmg; layout, kwargs...)

    if add_colorbars
        @show add_colorbars
        if p.edge_color[] != Makie.Automatic()
            cb = Colorbar(fig[1, 2], p.arrow_plot[])
        end
    end

    Makie.FigureAxisPlot(fig, ax, p)
end

function plot_nmg_acs(nmg; kwargs...)
    plot_nmg(nmg; get_ac_coloring_kwargs_simple(nmg)..., add_colorbars=false, kwargs...)
end

function make_full_reduced_singlestate_compgraph(L, codes;
    transprob_threshold=1000 * eps(),
    keep_max_only=false,
)
    nmgs = [make_ned_meta_graph(L, code; merge_shifted=true) for code in codes]

    # keys are those states which are a single state ac for some rule, values
    # are those rules which have it as a single-state ac
    ssac_states_to_codes = Dict{Vector{Int},Vector{Int}}()

    for ci in 1:length(codes)
        nmg = nmgs[ci]
        acs = attracting_components(nmg)
        for ac in acs
            # ignore non-single-state acs at this point
            if length(ac) == 1
                ssac_code = ac[1]
                ssac = label_for(nmg, ssac_code)

                c = codes[ci]
                if haskey(ssac_states_to_codes, ssac)
                    push!(ssac_states_to_codes[ssac], c)
                else
                    ssac_states_to_codes[ssac] = [c]
                end
            end
        end
    end
    ssac_states = sort(collect(keys(ssac_states_to_codes)))
    numssacs = length(ssac_states)

    cgraph = MetaGraph(
        DiGraph();
        label_type=Vector{Int},
        vertex_data_type=String,
        edge_data_type=Vector{Tuple{Int,Float64}},
        graph_data=(; L, codes, desc="reduced ma comp graph")
    )
    # add all the acs along with their labels
    for ssac in ssac_states
        cgraph[ssac] = join(ssac)
    end

    splitprobers = make_splitprober_mg.(nmgs)
    for src_ssac in ssac_states
        for (ci, c) in enumerate(codes)
            if c in ssac_states_to_codes[src_ssac]
                continue # skip self-edges
            end
            for dst_ssac in ssac_states
                if dst_ssac != src_ssac # skip self-edges
                    # now both are guaranteed to be acs of the given rule
                    # and we just need to find the probability of the system
                    # moving from src_ac to dst_ac under the given rule
                    tp = splitprobers[ci](src_ssac, dst_ssac)
                    if tp > transprob_threshold
                        if haskey(cgraph, src_ssac, dst_ssac)
                            push!(cgraph[src_ssac, dst_ssac], (c, tp))
                        else
                            cgraph[src_ssac, dst_ssac] = [(c, tp)]
                        end
                    end
                end
            end
        end
    end

    if !keep_max_only
        cgraph
    else
        reduce_allcodes_compgraph(cgraph)
    end
end

function condense_acs_mg(mg::MetaGraph{L,V,E,G}) where {L,V,E,G}
    all_acs = attracting_components(mg)
    nonsingle_ac_codes = filter(ac -> length(ac) != 1, all_acs)

    _taken_codes = reduce(vcat, nonsingle_ac_codes; init=Int[])

    # these will be the vertices in the partially condensed graph
    code_groups = []
    for v in 1:nv(mg)
        if !(v in _taken_codes)
            push!(code_groups, [v])
        end
    end
    code_groups = vcat(code_groups, nonsingle_ac_codes)
    label_groups = [label_for.(Ref(mg), cg) for cg in code_groups]

    # return label_groups

    sc_graph = SimpleWeightedDiGraph{Int,Float64}(length(code_groups))
    for (g1i, g1) in enumerate(label_groups)
        # can ignore non-single groups as they are acs and so have no outgoing edges
        if length(g1) == 1
            label1 = g1[1]
            for (g2i, g2) in enumerate(label_groups)
                if g1i == g2i
                    continue # tbf there should be none so this shoulnd't be needed
                end
                total_weight = 0.0
                for label2 in g2
                    if haskey(mg, label1, label2)
                        total_weight += mg[label1, label2]
                    end
                end
                if total_weight > 0.0
                    prev_weight = get_weight(sc_graph, g1i, g2i)
                    add_edge!(sc_graph, g1i, g2i, prev_weight + total_weight)
                end
            end
        end
    end

    label_groups, sc_graph
end

function make_splitprober_mg(g; error_if_not_ssac=false)
    label_groups, sc_graph = condense_acs_mg(g)
    group_is_ac = falses(length(label_groups))
    for gi in 1:nv(sc_graph)
        if isempty(outneighbors(sc_graph, gi))
            group_is_ac[gi] = true
        end
    end

    T = normalized_adjacency_matrix(sc_graph)
    lookup_matrix = inv(I(nv(sc_graph)) - T)

    function (l1, l2) # going from v1 to v2
        end_group_i = findfirst(x -> l2 in x, label_groups)
        if !group_is_ac[end_group_i] || (length(label_groups[end_group_i]) > 1)
            if error_if_not_ssac
                throw(ArgumentError("The end group is not a single-state ac"))
            else
                # 0. chance of splitting off to a non-ac state
                # and we disregard chances of reaching a particular state in a non-single-state ac as it's unclear
                return 0.0
            end
        end

        start_group_i = findfirst(x -> l1 in x, label_groups)
        if group_is_ac[start_group_i]
            return start_group_i == end_group_i ? 1.0 : 0.0
        else
            # so now we know the end state is a valid single-state ac
            # and that the start state is not a part of an ac
            return lookup_matrix[start_group_i, end_group_i]
        end
    end
end

function reduce_allcodes_compgraph_v2(
    accg::MetaGraph{C,G,L,V,Vector{Tuple{Int,Float64}}},
    threshold=1000 * eps()
) where {C,G,L,V}
    @show G
    accg_md = accg[]
    scgraph = MetaGraph(
        G();
        label_type=L,
        vertex_data_type=V,
        edge_data_type=Tuple{Vector{Int},Float64},
        graph_data=(; L=accg_md.L, codes=accg_md.codes, desc="reduced ma comp graph")
    )
    for l in labels(accg)
        scgraph[l] = accg[l]
    end
    for (s, d) in edge_labels(accg)
        full = accg[s, d]
        maxp = maximum(x -> x[2], full)
        goodenough_codes = []
        goodenough_vals = []
        for (c, v) in full
            if v >= maxp - threshold
                push!(goodenough_codes, c)
                push!(goodenough_vals, v)
            end
        end
        codes = sort(goodenough_codes)
        v = sum(goodenough_vals) / length(goodenough_vals)

        scgraph[s, d] = (codes, v)
    end
    scgraph
end

function gv_compgraph_sccs_v2(cg;
    elabel=:first,
    show_cond=true,
    filter_sccs=nothing,
    kwargs...
)
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        compound="true",
        kwargs...
    )

    code_sccs = strongly_connected_components(cg)
    cond = condensation(cg, code_sccs)
    sccs = [map(c -> label_for(cg, c), code_scc) for code_scc in code_sccs]

    if isnothing(filter_sccs)
        used_sscs = 1:length(sccs)
    else
        if isa(filter_sccs, Integer)
            xx = filter_sccs
            filter_sccs = x -> length(x) >= xx
        end
        used_sscs = findall(filter_sccs, sccs)
    end

    cluster_names = [(@sprintf "cluster %d" i) for i in 1:length(used_sscs)]
    clusters = [subgraph(g, cname) for cname in cluster_names]
    for (i, scc) in enumerate(sccs[used_sscs])
        for v in scc
            clusters[i] |> node(cg[v])
        end
    end
    fake_cluster_nodes = [(@sprintf "fcnode %d" ci) for ci in 1:length(clusters)]
    for (fcn, c) in zip(fake_cluster_nodes, clusters)
        c |> node(fcn;
            shape="point",
            style="invis"
        )
    end

    for scc in sccs[used_sscs]
        for v1 in scc
            for v2 in scc
                if v1 == v2
                    continue
                end
                if haskey(cg, v1, v2)
                    eval = cg[v1, v2]
                    if elabel == :first
                        label = string(eval[1][1])
                    elseif elabel == :all
                        label = join(string.(eval[1]), ",")
                    elseif elabel == :firstwithp
                        label = (@sprintf "%d, %.3g" eval[1][1] eval[2])
                    end
                    g |> edge(cg[v1], cg[v2]; label)
                end
            end
        end
    end

    if show_cond
        for (c1i, scc1i) in enumerate(used_sscs)
            for (c2i, scc2i) in enumerate(used_sscs)
                if scc1i == scc2i
                    continue
                end
                if has_edge(cond, scc1i, scc2i)
                    c1name = cluster_names[c1i]
                    c2name = cluster_names[c2i]
                    g |> edge(
                        fake_cluster_nodes[c1i],
                        fake_cluster_nodes[c2i];
                        ltail=c1name,
                        lhead=c2name,
                        # style="invis",
                        # weight="100"
                    )
                end
            end
        end
    end

    g
end

function gv_compgraph_v2(cg;
    elabel=:first,
    cluster_by=nothing,
    kwargs...
)
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    if isnothing(cluster_by)
        for v in labels(cg)
            g |> node(cg[v])
        end
    elseif cluster_by == :numones
        clusters = Dict{Int,Tuple{Any,Vector{Int}}}()
        for v in labels(cg)
            numones = count('1', cg[v])
            if !haskey(clusters, numones)
                clusters[numones] = (subgraph(g,
                        (@sprintf "cluster %d" numones),
                    ), [v])
            end
            xx = clusters[numones]
            xx[1] |> node(cg[v])
            push!(xx[2], v)
        end
        numspresent = sort(collect(keys(clusters)))
        for i in 1:(length(numspresent)-1)
            src_cluster_vs = clusters[numspresent[i]][2]
            dst_cluster_vs = clusters[numspresent[i+1]][2]
            for src_v in src_cluster_vs
                for dst_v in dst_cluster_vs
                    g |> edge(cg[src_v], cg[dst_v];
                        style="invis",
                        weight="100"
                    )
                end
            end
        end
    elseif cluster_by == :sccs
        sccs = strongly_connected_components(cg)
        clusters = [subgraph(g, (@sprintf "cluster %d" i)) for i in 1:length(sccs)]
        for (i, scc) in enumerate(sccs)
            for v in scc
                clusters[i] |> node(cg[v])
            end
        end
    else
        throw(ArgumentError("Unknown cluster_by value: $cluster_by"))
    end

    for (src, dst) in collect(edge_labels(cg))
        eval = cg[src, dst]
        if elabel == :first
            label = string(eval[1][1])
        elseif elabel == :all
            label = join(string.(eval[1]), ",")
        elseif elabel == :firstwithp
            label = (@sprintf "%d, %.3g" eval[1][1] eval[2])
        end
        g |> edge(cg[src], cg[dst]; label)
    end

    g
end


# FIX: Don't use, misses soem edges!
function merge_compgraph(cg)
    L = cg[].L

    itos(i) = NonEqDigits.itostate_safe(i, 2, L)

    mcg = MetaGraph(
        DiGraph();
        label_type=Int,
        vertex_data_type=String,
        edge_data_type=Tuple{Int,Float64},
        graph_data="ma comp graph"
    )

    unique_states = []
    unique_state_ids = []
    sti_to_usti = Dict{Int,Int}()
    for sti in labels(cg)
        st = itos(sti)
        usti = findfirst(ust -> states_related_by_shift(ust, st), unique_states)
        if isnothing(usti)
            push!(unique_state_ids, sti)
            push!(unique_states, st)
            usti = length(unique_state_ids)
        end
        sti_to_usti[sti] = usti
    end

    for usti in unique_state_ids
        mcg[usti] = cg[usti]
    end

    for (src_sti, dst_sti) in edge_labels(cg)
        src_usti = unique_state_ids[sti_to_usti[src_sti]]
        dst_usti = unique_state_ids[sti_to_usti[dst_sti]]

        if !haskey(mcg, src_usti, dst_usti)
            mcg[src_usti, dst_usti] = cg[src_sti, dst_sti]
        else
            (ne_val, ne_code) = cg[src_sti, dst_sti]
            prev = mcg[src_usti, dst_usti]
            if ne_val > prev[1]
                mcg[src_usti, dst_usti] = (ne_val, ne_code)
            end
        end
    end

    # unique_state_ids, sti_to_usti
    mcg
end

################################################################################
# Version 2
################################################################################
function make_full_singlestate_compgraph(L, codes;
    transprob_threshold=1000 * eps(),
    keep_max_only=true,
)
    numcodes = length(codes)
    neds = [make_ca_ned(L, c; show=false) for c in codes]
    graphs = graph.(neds)
    local numstates = GraphModels.numstates(neds[1])

    ss_acs_dict = Dict{Int,Vector{Int}}()
    for ci in 1:numcodes
        acs = attracting_components(graphs[ci])
        for ac in acs
            # ignore non-single-state acs at this point
            if length(ac) == 1
                ac = ac[1]

                c = codes[ci]
                if haskey(ss_acs_dict, ac)
                    push!(ss_acs_dict[ac], c)
                else
                    ss_acs_dict[ac] = [c]
                end
            end
        end
    end
    ss_acs = sort(collect(keys(ss_acs_dict)))
    numacs = length(ss_acs)

    cgraph = MetaGraph(
        DiGraph();
        label_type=Int,
        vertex_data_type=String,
        edge_data_type=Vector{Tuple{Int,Float64}},
        graph_data=(; L, codes, desc="ma comp graph")
    )
    # add all the acs along with their labels
    for ac in ss_acs
        cgraph[ac] = join(NonEqDigits.itostate_safe(ac, 2, L))
    end

    splitprobers = make_splitprober.(graphs)
    for src_ac in ss_acs
        for (ci, c) in enumerate(codes)
            if c in ss_acs_dict[src_ac]
                continue # skip self-edges
            end
            for dst_ac in ss_acs
                if dst_ac != src_ac # skip self-edges
                    # now both are guaranteed to be acs of the given rule
                    # and we just need to find the probability of the system
                    # moving from src_ac to dst_ac under the given rule
                    tp = splitprobers[ci](src_ac, dst_ac)
                    if tp > transprob_threshold
                        if haskey(cgraph, src_ac, dst_ac)
                            push!(cgraph[src_ac, dst_ac], (c, tp))
                        else
                            cgraph[src_ac, dst_ac] = [(c, tp)]
                        end
                    end
                end
            end
        end
    end

    if !keep_max_only
        cgraph
    else
        reduce_allcodes_compgraph_v2(cgraph)
    end
end

function reduce_allcodes_compgraph(accg::MetaGraph{C,G,L,V,Vector{Tuple{Int,Float64}}}) where {C,G,L,V}
    accg_md = accg[]
    scgraph = MetaGraph(
        G();
        label_type=L,
        vertex_data_type=V,
        edge_data_type=Tuple{Int,Float64},
        graph_data=(; L=accg_md.L, codes=accg_md.codes, desc="reduced ma comp graph")
    )
    for l in labels(accg)
        scgraph[l] = accg[l]
    end
    for (s, d) in edge_labels(accg)
        full = accg[s, d]
        maxi = findmax(x -> x[2], full)[2]
        scgraph[s, d] = full[maxi]
    end
    scgraph
end

# Figuring out splitting probabilities
function condense_acs(g::SimpleWeightedDiGraph)
    all_acs = attracting_components(g)
    nonsingle_acs = filter(ac -> length(ac) != 1, all_acs)

    _taken_vs = reduce(vcat, nonsingle_acs; init=Int[])

    # these will be the vertices in the partially condensed graph
    groups = []
    for v in 1:nv(g)
        if !(v in _taken_vs)
            push!(groups, [v])
        end
    end
    groups = vcat(groups, nonsingle_acs)

    v_to_group = fill(-1, nv(g))
    for v in 1:nv(g)
        v_to_group[v] = findfirst(x -> v in x, groups)
    end

    sc_graph = SimpleWeightedDiGraph{Int,Float64}(length(groups))
    numgroups = length(groups)
    for (g1i, g1) in enumerate(groups)
        # can ignore non-single groups as they are acs and so have no outgoing edges
        if length(g1) == 1
            v1 = g1[1]
            for (g2i, g2) in enumerate(groups)
                if g1i == g2i
                    continue # tbf there should be none so this shoulnd't be needed
                end
                total_weight = 0.0
                for v2 in g2
                    total_weight += get_weight(g, v1, v2)
                end
                if total_weight > 0.0
                    prev_weight = get_weight(sc_graph, g1i, g2i)
                    add_edge!(sc_graph, g1i, g2i, prev_weight + total_weight)
                end
            end
        end
    end

    groups, sc_graph
end

function viz_condense_acs(g)
    groups, sc_graph = condense_acs(g)
    f = Figure()
    layout = Spring(dim=3)
    ax1 = LScene(f[1, 1])
    ax2 = LScene(f[1, 2])
    graphplot!(ax1, g;
        nlabels=string.(1:nv(g)),
        elabels=[(@sprintf "%.3g" weight(e)) for e in edges(g)],
        layout,
        get_ac_coloring_kwargs(g)...
    )
    graphplot!(ax2, sc_graph;
        nlabels=repr.(groups),
        elabels=[(@sprintf "%.3g" weight(e)) for e in edges(sc_graph)],
        layout
    )
    FigureAxisAnything(f, [ax1, ax2], (; groups, sc_graph))
end

"""
Returns an anonymous function that can be passed two states of the graph g
the second of which has to be a single-state ac. Will return the probability
that if started in the first state that the system will end up in the second state.
"""
function make_splitprober(g; error_if_not_ssac=false)
    groups, sc_graph = condense_acs(g)
    group_is_ac = falses(length(groups))
    for gi in 1:nv(sc_graph)
        if isempty(outneighbors(sc_graph, gi))
            group_is_ac[gi] = true
        end
    end

    T = normalized_adjacency_matrix(sc_graph)
    lookup_matrix = inv(I(nv(sc_graph)) - T)

    function (v1, v2) # going from v1 to v2
        start_group_i = findfirst(x -> v1 in x, groups)
        end_group_i = findfirst(x -> v2 in x, groups)

        if !group_is_ac[end_group_i] || (length(groups[end_group_i]) > 1)
            if error_if_not_ssac
                throw(ArgumentError("The end group is not a single-state ac"))
            else
                return 0.0
            end
        elseif group_is_ac[start_group_i]
            return start_group_i == end_group_i ? 1.0 : 0.0
        else
            # so now we know the end state is a valid single-state ac
            # and that the start state is not a part of an ac
            return lookup_matrix[start_group_i, end_group_i]
        end
    end
end

function make_splitprober_full(g)
    groups, sc_graph = condense_acs(g)
    group_is_ac = falses(length(groups))
    for gi in 1:nv(sc_graph)
        if isempty(outneighbors(sc_graph, gi))
            group_is_ac[gi] = true
        end
    end

    T = normalized_adjacency_matrix(sc_graph)
    lookup_matrix = inv(I(nv(sc_graph)) - T)

    function (v) # v is a vertex of the g graph in which the system starts
        start_group = findfirst(x -> v in x, groups)
        if group_is_ac[start_group]
            return [(1.0, groups[start_group])]
        end
        group_probs = @view lookup_matrix[start_group, :]
        results = Tuple{Float64,Vector{Int}}[]
        # for (group, group_prob) in zip(groups, group_probs)
        for gi in 1:length(groups)
            if group_is_ac[gi]
                group = groups[gi]
                group_prob = group_probs[gi]
                if group_prob > 0.0
                    push!(results, (group_prob, group))
                end
            end
        end
        return results
    end
end

# v2 plotting
function plot_compgraph(cg;
    highlight_reachable=true,
    highlight_reversibility=true,
    base_color=:black,
    include_state_code=true,
    elabels=nothing,
    kwargs...
)
    auto_kwargs = Dict{Symbol,Any}()
    auto_kwargs[:node_size] = 15

    if highlight_reachable == true
        highlight_reachable = :orangered3
    end
    if !isnothing(highlight_reachable)
        auto_kwargs[:node_color] = [isempty(inneighbors(cg, v)) ? base_color : highlight_reachable for v in 1:nv(cg)]
    end

    if highlight_reversibility
        colors = []
        for (src, dst) in edge_labels(cg)
            color = base_color
            if haskey(cg, dst, src)
                color = ColorSchemes.managua10[5]
            else
                src_numones = count('1', cg[src])
                dst_numones = count('1', cg[dst])
                if src_numones > dst_numones
                    color = ColorSchemes.managua10[1]
                else
                    color = ColorSchemes.managua10[10]
                end
            end

            push!(colors, color)
        end
        auto_kwargs[:edge_color] = colors
    end

    auto_kwargs[:nlabels] = if include_state_code
        [(@sprintf "%s(%d)" cg[l] l) for l in labels(cg)]
    else
        string.(getindex.(Ref(cg), labels(cg)))
    end

    if elabels == :codeonly
        auto_kwargs[:elabels] = [string(cg[s, d][1]) for (s, d) in edge_labels(cg)]
    end

    graphplot(cg;
        auto_kwargs..., kwargs...
    )
end

function highlight_sccs_kwargs(g; highlight_color=:red, base_color=:black)
    sccs = filter(x -> length(x) > 1, strongly_connected_components(g))

    node_color = [base_color for _ in 1:nv(g)]
    for scc in sccs
        for v in scc
            node_color[v] = highlight_color
        end
    end

    edge_color = [base_color for _ in 1:ne(g)]
    for (ei, e) in enumerate(edges(g))
        src_scc = findfirst(x -> e.src in x, sccs)
        dst_scc = findfirst(x -> e.dst in x, sccs)
        if !isnothing(src_scc) && !isnothing(dst_scc)
            if src_scc == dst_scc
                edge_color[ei] = highlight_color
            end
        end
    end

    (; node_color, edge_color)
end

function gv_compgraph(cg;
    cluster_by=nothing,
    kwargs...
)
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    if isnothing(cluster_by)
        for v in labels(cg)
            g |> node(cg[v])
        end
    elseif cluster_by == :numones
        clusters = Dict{Int,Tuple{Any,Vector{Int}}}()
        for v in labels(cg)
            numones = count('1', cg[v])
            if !haskey(clusters, numones)
                clusters[numones] = (subgraph(g,
                        (@sprintf "cluster %d" numones),
                    ), [v])
            end
            xx = clusters[numones]
            xx[1] |> node(cg[v])
            push!(xx[2], v)
        end
        numspresent = sort(collect(keys(clusters)))
        for i in 1:(length(numspresent)-1)
            src_cluster_vs = clusters[numspresent[i]][2]
            dst_cluster_vs = clusters[numspresent[i+1]][2]
            for src_v in src_cluster_vs
                for dst_v in dst_cluster_vs
                    g |> edge(cg[src_v], cg[dst_v];
                        style="invis",
                        weight="100"
                    )
                end
            end
        end
    elseif cluster_by == :sccs
        sccs = strongly_connected_components(cg)
        clusters = [subgraph(g, (@sprintf "cluster %d" i)) for i in 1:length(sccs)]
        for (i, scc) in enumerate(sccs)
            for v in scc
                clusters[i] |> node(cg[v])
            end
        end
    else
        throw(ArgumentError("Unknown cluster_by value: $cluster_by"))
    end

    for (src, dst) in collect(edge_labels(cg))
        g |> edge(cg[src], cg[dst];
            label=(@sprintf "%d" cg[src, dst][1]),
            # weight="10",
        )
    end

    g
end

function gv_compgraph_sccs(cg;
    show_cond=true,
    skip_singles=false,
    kwargs...
)
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        compound="true",
        kwargs...
    )

    code_sccs = strongly_connected_components(cg)
    cond = condensation(cg, code_sccs)
    sccs = [map(c -> label_for(cg, c), code_scc) for code_scc in code_sccs]

    used_sscs = 1:length(sccs)
    if skip_singles
        used_sscs = findall(x -> length(x) > 1, sccs)
    end

    cluster_names = [(@sprintf "cluster %d" i) for i in 1:length(used_sscs)]
    clusters = [subgraph(g, cname) for cname in cluster_names]
    for (i, scc) in enumerate(sccs[used_sscs])
        for v in scc
            clusters[i] |> node(cg[v])
        end
    end
    fake_cluster_nodes = [(@sprintf "fcnode %d" ci) for ci in 1:length(clusters)]
    for (fcn, c) in zip(fake_cluster_nodes, clusters)
        c |> node(fcn;
            shape="point",
            style="invis"
        )
    end

    for scc in sccs[used_sscs]
        for v1 in scc
            for v2 in scc
                if v1 == v2
                    continue
                end
                if haskey(cg, v1, v2)
                    g |> edge(cg[v1], cg[v2];
                        label=(@sprintf "%d" cg[v1, v2][1]),
                        # weight="10",
                    )
                end
            end
        end
    end

    if show_cond
        for (c1i, scc1i) in enumerate(used_sscs)
            for (c2i, scc2i) in enumerate(used_sscs)
                if scc1i == scc2i
                    continue
                end
                if has_edge(cond, scc1i, scc2i)
                    c1name = cluster_names[c1i]
                    c2name = cluster_names[c2i]
                    g |> edge(
                        fake_cluster_nodes[c1i],
                        fake_cluster_nodes[c2i];
                        ltail=c1name,
                        lhead=c2name,
                        # style="invis",
                        # weight="100"
                    )
                end
            end
        end
    end

    g
end

################################################################################
# Version 1, don't use
# Makes single-state ac graphs for each rule separately
################################################################################
function make_singlestate_ac_compgraphs(L, codes;
    prob_threshold=1000 * eps()
)
    numcodes = length(codes)
    neds = [make_ca_ned(L, c; show=false) for c in codes]
    graphs = graph.(neds)
    local numstates = GraphModels.numstates(neds[1])

    delta = I(numstates)

    nams = normalized_adjacency_matrix.(neds)
    ac_to_ac_pmats = [inv(delta - nam) for nam in nams]

    allacs = Dict{Int,Vector{Int}}()
    for ci in 1:numcodes
        acs = attracting_components(graphs[ci])
        for ac in acs
            if length(ac) != 1
                throw(ArgumentError("Found a non single-state ac"))
            end
            ac = ac[1]

            if haskey(allacs, ac)
                push!(allacs[ac], ci)
            else
                allacs[ac] = [ci]
            end
        end
    end
    numacs = length(allacs)

    # make a separate graph for each rule, then potentially merge them later
    sccg_template = MetaGraph(
        DiGraph();
        label_type=Int,
        vertex_data_type=String,
        edge_data_type=Float64,
        graph_data="single rule ma comp graph"
    )
    for (ac, _) in allacs
        sccg_template[ac] = join(NonEqDigits.itostate_safe(ac, 2, L))
    end
    sccgs = [copy(sccg_template) for _ in 1:numcodes]

    for (ci, c) in enumerate(codes)
        cacg = sccgs[ci]
        for (ac, accodeis) in allacs
            if ci in accodeis # skip those with this ac (they would be self-edges)
                continue
            end
            # @show ac
            for dst_ac in keys(allacs)
                if dst_ac == ac
                    continue
                end
                trans_prob = ac_to_ac_pmats[ci][ac, dst_ac]
                if abs(trans_prob) > prob_threshold
                    @show trans_prob
                    add_edge!(cacg, ac, dst_ac, trans_prob)
                end
            end
        end
    end

    keys(allacs), sccgs
end

function merge_sccgs(L, codes, sccgs; acs=collect(labels(sccgs[1])))
    xx = MetaGraph(
        DiGraph();
        label_type=Int,
        edge_data_type=Dict{Int,Float64},
        graph_data="temp meta graph"
        # graph_data="merged (many rule) ma comp graph"
    )
    # make sure we have all the acs in the graph
    for ac in acs
        xx[ac] = nothing
    end

    for (c, sccg) in zip(codes, sccgs)
        for (src_ac, dst_ac) in edge_labels(sccg)
            trans_prob = sccg[src_ac, dst_ac]
            if !haskey(xx, src_ac, dst_ac)
                xx[src_ac, dst_ac] = Dict{Int,Float64}(c => trans_prob)
            else
                xx[src_ac, dst_ac][c] = trans_prob
            end
        end
    end

    yy = MetaGraph(
        DiGraph();
        label_type=Int,
        vertex_data_type=String,
        edge_data_type=Tuple{Int,Float64},
        graph_data="temp meta graph"
        # graph_data="merged (many rule) ma comp graph"
    )
    for ac in acs
        yy[ac] = join(NonEqDigits.itostate_safe(ac, 2, L))
    end
    for (src_ac, dst_ac) in edge_labels(xx)
        dd = xx[src_ac, dst_ac]
        max_trans_prob, max_code = findmax(dd)
        yy[src_ac, dst_ac] = (max_code, max_trans_prob)
    end

    yy
end

function single_rule_comp_gv(L, acs, acg;
    force_order=true,
    include_isolated=true,
    double_arrows=false,
    kwargs...
)
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    node_clusters = [subgraph(g, (@sprintf "cluster %s" string(i))) for i in 0:L]
    cluster_is = [[] for _ in 0:L]
    for ac in acs
        st = NonEqDigits.itostate_safe(ac, 2, L)
        num_ones = count(x -> x == 1, st)
        cluster_index = 1 + num_ones
        node_clusters[cluster_index] |> node(string(ac); label=join(st))
        push!(cluster_is[cluster_index], ac)
    end

    nonempty_cluster_is = filter(!isempty, cluster_is)
    if force_order
        for ci in 1:(L-1)
            for i1 in nonempty_cluster_is[ci]
                for i2 in nonempty_cluster_is[ci+1]
                    g |> edge(string(i1), string(i2);
                        style="invis",
                        # weight="100"
                    )
                end
            end
        end
    end

    if double_arrows
        ng = graph(ned)
        for i in 1:numstates(ned)
            for j in (i+1):numstates(ned)
                itoj = !iszero(get_weight(ng, i, j))
                jtoi = !iszero(get_weight(ng, j, i))
                if itoj && jtoi
                    g |> edge(string(i), string(j);
                        dir="both",
                        # weight="10"
                    )
                elseif itoj
                    g |> edge(
                        string(i),
                        string(j);
                        # weight="10"
                    )
                elseif jtoi
                    g |> edge(
                        string(j),
                        string(i);
                        # weight="10"
                    )
                end
            end
        end
    else
        for e in edges(acg)
            src = acs[e.src]
            dst = acs[e.dst]
            g |> edge(
                string(src),
                string(dst);
                # weight="10"
            )
        end
    end

    g
end
