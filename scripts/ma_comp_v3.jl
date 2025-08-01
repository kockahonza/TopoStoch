using DrWatson;
@quickactivate "TopoStochSim";

using GraphModels
using NonEqDigits
using ColorSchemes
using GraphvizDotLang: digraph, subgraph, node, edge
using MetaGraphsNext
using DataFrames

using TernaryDiagrams
using PythonCall

################################################################################
# Util
################################################################################
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

# general metagraph plotting
function metagraph_kwargs(mg::MetaGraph{C,G,L,VD,ED}) where {C,G,L,VD,ED}
    rslt = (;
        nlabels=string.(getindex.(Ref(mg), labels(mg))),
    )
    if ED <: Number
        rslt = (;
            elabels=[(@sprintf "%.3g" mg[src, dst]) for (src, dst) in edge_labels(mg)],
            edge_color=[mg[src, dst] for (src, dst) in edge_labels(mg)],
            rslt...
        )
    else
        rslt = (;
            elabels=[string(mg[src, dst]) for (src, dst) in edge_labels(mg)],
            rslt...
        )
    end
    rslt
end
function metagraphplot(mg; kwargs...)
    graphplot(mg; metagraph_kwargs(mg)..., kwargs...)
end
function metagraphplot!(ax, mg; kwargs...)
    graphplot!(ax, mg; metagraph_kwargs(mg)..., kwargs...)
end

################################################################################
# Making a transition graph for a given rule
################################################################################
"""
This makes a graph of all states with transitions between them that result from
the given rule. The result is a metagraph with the states (Int vectors) as vertex
labels and floats as edge values/weights. Stores the code and L in metadata.
If `merge_shifted` is set to true, then the states that are related by shifts
are merged into a single state.
"""
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

# Some plotting for these
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
        if p.edge_color[] isa Vector
            cb = Colorbar(fig[1, 2], p.arrow_plot[])
        else
            @show typeof(p.edge_color[])
            @debug "skipping colorbars as edge_color is a $(typeof(p.edge_color[]))"
        end
    end

    Makie.FigureAxisPlot(fig, ax, p)
end

function plot_nmg_acs(nmg; kwargs...)
    plot_nmg(nmg; get_ac_coloring_kwargs_simple(nmg)..., add_colorbars=false, kwargs...)
end

################################################################################
# Making the computational graphs
################################################################################
function make_ssac_compgraphs(L, codes;
    transprob_threshold=1000 * eps(),
    merge_shifted=false,
    include_self_edges=true
)
    nmgs = [make_ned_meta_graph(L, code; merge_shifted) for code in codes]

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

    if isa(transprob_threshold, Number)
        _xax = transprob_threshold
        transprob_threshold = x -> x > _xax
    end

    splitprobers = make_ssac_splitprober.(nmgs)
    for src_ssac in ssac_states
        for (ci, c) in enumerate(codes)
            if c in ssac_states_to_codes[src_ssac]
                if include_self_edges
                    if haskey(cgraph, src_ssac, src_ssac)
                        push!(cgraph[src_ssac, src_ssac], (c, 1.0))
                    else
                        cgraph[src_ssac, src_ssac] = [(c, 1.0)]
                    end
                end
                continue
            end
            for dst_ssac in ssac_states
                if dst_ssac != src_ssac # skip self-edges
                    # now both are guaranteed to be acs of the given rule
                    # and we just need to find the probability of the system
                    # moving from src_ac to dst_ac under the given rule
                    tp = splitprobers[ci](src_ssac, dst_ssac)
                    if transprob_threshold(tp)
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

    cgraph
end

"""
Returns an anonymous function that can be passed two states of the graph g
the second of which has to be a single-state ac. Will return the probability
that if started in the first state that the system will end up in the second state.
"""
function make_ssac_splitprober(g; error_if_not_ssac=false)
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

################################################################################
# Marge all single rule compgraphs into one, keeping a list of rules per transition
################################################################################
function merge_ssac_compgraphs(
    accg::MetaGraph{C,G,L,V,Vector{Tuple{Int,Float64}}},
    threshold=1000 * eps()
) where {C,G,L,V}
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

function make_merged_ssac_compgraph(args...; kwargs...)
    merge_ssac_compgraphs(make_ssac_compgraphs(args...; kwargs...))
end

################################################################################
# Analysis on compgraphs
################################################################################
"""Works on a merged compgraph"""
function find_rule_groups(mcg)
    groups = []
    for (s, d) in edge_labels(mcg)
        e = mcg[s, d]
        erules = sort(e[1])
        if !(erules in groups)
            push!(groups, erules)
        end
    end
    sort!(groups)
    sort!(groups; by=length)
end

function find_outcomes(mcg, start_label, rule)
    results = []
    for dst in outneighbor_labels(mcg, start_label)
        e = mcg[start_label, dst]
        if rule in e[1]
            push!(results, (dst, e[2]))
        end
    end
    results
end
function find_unique_outcome(mcg, start_label, rule)
    outcomes = find_outcomes(mcg, start_label, rule)
    if isempty(outcomes)
        nothing
    elseif length(outcomes) == 1
        if abs(outcomes[1][2] - 1.0) > 1000 * eps()
            @warn "getting unique outcomes that are not with 100% probablity"
        end
        outcomes[1][1]
    else
        throw(ArgumentError("Multiple outcomes found for rule $rule from state $start_label: $outcomes"))
    end
end

"""Works on a merged compgraph"""
function test_commutative(mcg)
    rule_groups = find_rule_groups(mcg)
    present_rules = sort(unique(reduce(vcat, rule_groups)))
    all_states = labels(mcg)

    # Possible outcomes
    # - they both lead to nothing
    # - one leads to nothing
    # - both non-nothing and equal
    # - both non-nothing and different

    full_rslt = DataFrame(;
        r1=Int[],
        r2=Int[],
        st1=String[],
        code=Int[],
        f12=Union{Nothing,String}[],
        f21=Union{Nothing,String}[],
        final_state=Union{Nothing,String}[],
    )

    ronly_rslts = DataFrame(;
        r1=Int[],
        r2=Int[],
        alwayscommutative=Bool[],
        dosomething=Bool[],
    )

    for r1 in present_rules
        for r2 in present_rules
            r1r2always_commutative = true
            r1r2dosomething = false
            for st1 in all_states
                r1st2 = find_unique_outcome(mcg, st1, r1)
                r1final = isnothing(r1st2) ? nothing : find_unique_outcome(mcg, r1st2, r2)
                r2st2 = find_unique_outcome(mcg, st1, r2)
                r2final = isnothing(r2st2) ? nothing : find_unique_outcome(mcg, r2st2, r1)

                code = -1
                final_state = nothing
                if isnothing(r1final) && isnothing(r2final)
                    code = 1
                elseif isnothing(r1final) || isnothing(r2final)
                    r1r2always_commutative = false
                    r1r2dosomething = true
                    code = 2
                else
                    r1r2dosomething = true
                    if r1final != r2final
                        r1r2always_commutative = false
                        code = 3
                    else
                        code = 4
                        final_state = r1final
                    end
                end
                if !isnothing(final_state)
                    final_state = mcg[final_state]
                end
                if !isnothing(r1final)
                    r1final = mcg[r1final]
                end
                if !isnothing(r2final)
                    r2final = mcg[r2final]
                end
                push!(full_rslt, (r1, r2, mcg[st1], code, r1final, r2final, final_state))
            end
            push!(ronly_rslts, (r1, r2, r1r2always_commutative, r1r2dosomething))
        end
    end

    full_rslt, ronly_rslts
end

function find_noncommutative(mcg)
    rule_groups = find_rule_groups(mcg)
    present_rules = sort(unique(reduce(vcat, rule_groups)))
    all_states = labels(mcg)

    full_rslt = DataFrame(;
        r1=Int[],
        r2=Int[],
        src=String[],
        code=Int[],
        s1=Union{Nothing,String}[],
        s2=Union{Nothing,String}[],
        s12=Union{Nothing,String}[],
        s21=Union{Nothing,String}[],
    )

    for r1 in present_rules
        for r2 in present_rules
            for st1 in all_states
                s1 = find_unique_outcome(mcg, st1, r1)
                if isnothing(s1)
                    s12 = nothing
                else
                    s12 = find_unique_outcome(mcg, s1, r2)
                    s1 = mcg[s1]
                    if !isnothing(s12)
                        s12 = mcg[s12]
                    end
                end

                s2 = find_unique_outcome(mcg, st1, r2)
                if isnothing(s2)
                    s21 = nothing
                else
                    s21 = find_unique_outcome(mcg, s2, r1)
                    s2 = mcg[s2]
                    if !isnothing(s21)
                        s21 = mcg[s21]
                    end
                end

                code = -1
                if isnothing(s12) && isnothing(s21)
                    code = 1
                elseif isnothing(s12) || isnothing(s21)
                    code = 2
                else
                    if s12 != s21
                        code = 3
                    else
                        code = 4
                    end
                end
                push!(full_rslt, (r1, r2, mcg[st1], code, s1, s2, s12, s21))
            end
        end
    end

    full_rslt
end

################################################################################
# Vizualisation
################################################################################
function gv_compgraph_sccs(cg;
    elabel=:first,
    show_cond=true,
    filter_sccs=nothing,
    skip_self_edges=true,
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
                if skip_self_edges && (v1 == v2)
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
                    )
                end
            end
        end
    end

    g
end

function gv_compgraph(cg;
    elabel=:first,
    cluster_by=nothing,
    skip_self_edges=true,
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
        clusters = Dict{Int,Tuple{Any,Vector{Vector{Int}}}}()
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
                clusters[i] |> node(cg[label_for(cg, v)])
            end
        end
    else
        throw(ArgumentError("Unknown cluster_by value: $cluster_by"))
    end

    for (src, dst) in collect(edge_labels(cg))
        if skip_self_edges && (src == dst)
            continue
        end
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

function plot_compgraph(cg;
    highlight_reachable=true,
    highlight_reversibility=true,
    base_color=:black,
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

    auto_kwargs[:nlabels] = string.(getindex.(Ref(cg), labels(cg)))

    if elabels == :codeonly
        auto_kwargs[:elabels] = [string(cg[s, d][1]) for (s, d) in edge_labels(cg)]
    end

    graphplot(cg;
        auto_kwargs..., kwargs...
    )
end

################################################################################
# Basins of attractions general
################################################################################
function all_rules_with_all0sand1s(L, rules=0:255;
    onlyacs=false,
    merge_shifted=false,
    return_all=false
)
    all0s = fill(0, L)
    all1s = fill(1, L)
    rslt = []

    for r in rules
        nmg = make_ned_meta_graph(L, r; merge_shifted)
        if isempty(outneighbor_labels(nmg, all0s)) &&
           isempty(outneighbor_labels(nmg, all1s))
            if onlyacs
                acs = attracting_components(nmg)
                if length(acs) != 2
                    continue
                end
            end
            push!(rslt, r)
        end
    end

    ucodes = ca_ucodes_f1()
    sb = []
    for r in rslt
        eqp = find_eq_parent_ucode(r; ucodes)
        rne = ca_numenzymes(r)
        eqpne = ca_numenzymes(eqp)
        push!(sb, (eqpne, eqp, rne))
    end
    ii = sortperm(sb)

    if return_all
        rslt[ii], sb[ii]
    else
        rslt[ii]
    end
end

function get_all_ac_split_probs(nmg::MetaGraph{C,G,L}) where {C,G,L}
    ff = make_full_splitprober(nmg)
    d = Dict{L,Dict{Vector{L},Float64}}()
    for l in labels(nmg)
        d[l] = ff(l)
    end
    d
end

function make_full_splitprober(g; threshold=1000 * eps())
    label_groups, sc_graph = condense_acs_mg(g)
    group_is_ac = falses(length(label_groups))
    for gi in 1:nv(sc_graph)
        if isempty(outneighbors(sc_graph, gi))
            group_is_ac[gi] = true
        end
    end

    T = normalized_adjacency_matrix(sc_graph)
    lookup_matrix = inv(I(nv(sc_graph)) - T)

    function (st) # going from v1 to v2
        start_group_i = findfirst(x -> st in x, label_groups)
        if isnothing(start_group_i)
            throw(ArgumentError("The start state $st is not present"))
        end
        # @show start_group_i st label_groups
        if group_is_ac[start_group_i]
            Dict{eltype(label_groups),Float64}(label_groups[start_group_i] => 1.0)
        else
            rslt = Dict{eltype(label_groups),Float64}()
            for eg_i in 1:length(label_groups)
                if group_is_ac[eg_i]
                    p = lookup_matrix[start_group_i, eg_i]
                    if p > threshold
                        rslt[label_groups[eg_i]] = p
                    end
                end
            end
            total = sum(values(rslt))
            if abs(total - 1) > threshold
                @warn "Getting a full split prob vector with probability not equal to 1! The total is $total"
            end

            for k in keys(rslt)
                rslt[k] /= total
            end

            rslt
        end
    end
end

################################################################################
# All 0s vs all 1s, mostly plotting
################################################################################
function get_all0sv1s_colors(nmg::MetaGraph{C,G,LT};
    colormap=ColorSchemes.RdBu_4,
    colorrange=(0.0, 1.0),
    # colorscale=identity,
    checks=true,
) where {C,G,LT}
    if checks
        acs = attracting_components(nmg)
        if length(acs) != 2
            throw(ArgumentError("The given nmg does not have exactly two acs, it has $(length(acs))"))
        end
        hasall0s = false
        hasall1s = false
        for ac in acs
            if length(ac) != 1
                throw(ArgumentError("The given nmg does not have single-state acs, it has $ac"))
            end
            st = label_for(nmg, ac[1])
            if all(x -> x == 0, st)
                hasall0s = true
            end
            if all(x -> x == 1, st)
                hasall1s = true
            end
        end
        if !hasall0s || !hasall1s
            throw(ArgumentError("The given nmg does not have all 0s and all 1s states as acs"))
        end
    end

    L = nmg[].L
    all0sstate = fill(0, L)
    all1sstate = fill(1, L)

    color_vals = Dict{LT,Float64}()
    ff = make_full_splitprober(nmg)

    for l in labels(nmg)
        sps = ff(l)
        color_vals[l] = get(sps, [all0sstate], 0.0)

        if checks
            if length(keys(sps)) == 1
                ac = collect(keys(sps))[1]
                if (ac != [all0sstate]) && (ac != [all1sstate])
                    @warn "The state $l has split probabilities to $ac, which is not all 0s or all 1s"
                end
            elseif length(keys(sps)) == 2
                ac1 = collect(keys(sps))[1]
                ac2 = collect(keys(sps))[2]

                if !(
                    ((ac1 == [all0sstate]) && (ac2 == [all1sstate])) ||
                    ((ac1 == [all1sstate]) && (ac2 == [all0sstate]))
                )
                    @warn "The state $l has split probabilities to $ac1 and $ac2, which are not all 0s and all 1s"
                end
            else
                @warn "The state $l does not have exactly two split probabilities, it has $(length(keys(sps)))"
            end
        end
    end

    if isnothing(colorrange)
        mincv, maxcv = extrema(values(color_vals))
    else
        mincv, maxcv = colorrange
    end
    deltacv = maxcv

    colors = Dict{LT,Color}()
    for k in keys(color_vals)
        v = (color_vals[k] - mincv) / deltacv
        colors[k] = get(colormap, v)
    end

    colors
end

function get_all0sv1s_colors_v2(nmg::MetaGraph{C,G,LT};
    cmap0sv1s=ColorSchemes.RdBu_4,
    nondetcolor=ColorSchemes.viridis[end],
) where {C,G,LT}
    L = nmg[].L
    all0sstate = fill(0, L)
    all1sstate = fill(1, L)

    colors = Dict{LT,Color}()
    ff = make_full_splitprober(nmg)

    for l in labels(nmg)
        sps = ff(l)
        prob_all0s = get(sps, [all0sstate], 0.0)
        prob_all1s = get(sps, [all1sstate], 0.0)
        prob_nonextreme = 0.0
        for ac in keys(sps)
            if ac != [all0sstate] && ac != [all1sstate]
                prob_nonextreme += sps[ac]
            end
        end

        val_0sv1s = if prob_all0s == 0.0
            0.0
        else
            prob_all0s / (prob_all0s + prob_all1s)
        end

        # val_0sv1s = prob_all1s

        c0v1s = get(cmap0sv1s, val_0sv1s)
        cmap = ColorScheme(range(c0v1s, nondetcolor, 100))
        colors[l] = get(cmap, prob_nonextreme)
    end

    colors
end

function all0sv1s_colors_v2_colorbar(N=100;
    cmap0sv1s=ColorSchemes.RdBu_4,
    nondetcolor=ColorSchemes.viridis[end],
)
    probs_all0s = range(0.0, 1.0, N)
    probs_nonextreme = range(0.0, 1.0, N)

    function func(pall0s, pnonextreme)
        c0v1s = get(cmap0sv1s, pall0s)
        cmap = ColorScheme(range(c0v1s, nondetcolor, 100))
        get(cmap, pnonextreme)
    end

    fap = heatmap(probs_all0s, probs_nonextreme, func)
    fap.axis.xlabel = "probability of all 1s vs all 0s"
    fap.axis.ylabel = "probability of reaching any other ac"

    fap
end

function all0sv1s_colors_v2_colorbar_v2!(ax, N=100;
    epsilon=0.0,
    cmap0sv1s=ColorSchemes.RdBu_4,
    nondetcolor=ColorSchemes.viridis[end],
)
    function func(prob_all0s, prob_nonextreme)
        prob_all1s = 1.0 - prob_all0s - prob_nonextreme

        val_0sv1s = if prob_all0s == 0.0
            0.0
        else
            prob_all0s / (prob_all0s + prob_all1s)
        end

        # val_0sv1s = prob_all1s

        c0v1s = get(cmap0sv1s, val_0sv1s)
        cmap = ColorScheme(range(c0v1s, nondetcolor, 100))
        get(cmap, prob_nonextreme)
    end

    p0s = []
    p1s = []
    pes = []
    cval = []
    color = []

    for p0 in range(epsilon, 1.0 - epsilon, N)
        for p1 in range(epsilon, 1.0 - p0 - epsilon, N)
            pelse = 1.0 - p0 - p1
            if (pelse <= epsilon) || (pelse >= (1.0 - epsilon))
                continue
            end
            push!(p0s, p0)
            push!(p1s, p1)
            push!(pes, pelse)

            push!(cval, rand() * p1 - p0)
            push!(color, func(p0, pelse))
        end
    end

    ternaryscatter!(ax, p0s, p1s, pes; color=color, marker=:circle)
    # ternarycontourf!(ax, p0s, p1s, pes, cval; levels=10)

    ternaryaxis!(ax;
        labelx="p0",
        labely="p1",
        labelz="pelse",
    )
    xlims!(ax, -0.2, 1.2)
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax)

    nothing
end

function all0sv1s_colors_v2_colorbar_v5!(ax, N=20;
    cmap0sv1s=ColorSchemes.RdBu_4,
    nondetcolor=ColorSchemes.viridis[end],
)
    function func(prob_all0s, prob_nonextreme)
        prob_all1s = 1.0 - prob_all0s - prob_nonextreme

        val_0sv1s = if prob_all0s == 0.0
            0.0
        else
            prob_all0s / (prob_all0s + prob_all1s)
        end

        # val_0sv1s = prob_all1s

        c0v1s = get(cmap0sv1s, val_0sv1s)
        cmap = ColorScheme(range(c0v1s, nondetcolor, 100))
        get(cmap, prob_nonextreme)
    end

    # make the mesh
    vertices = []
    # centers = []
    colors = []
    vi_rows = []

    a = 1 / (N - 1)
    vi = 1
    for ri in 1:N
        push!(vi_rows, [])
        for j in 1:(N-ri+1)
            push!(vi_rows[ri], vi)

            y = (ri - 1) * (sqrt(3) / 2) * a
            x = (ri - 1) * (a / 2) + a * (j - 1)
            pos = [x, y]

            push!(vertices, pos)

            pe = (2 / sqrt(3)) * y
            p0 = 1 - x - y / sqrt(3)
            p1 = x - y / sqrt(3)


            # color = get(ColorSchemes.viridis, p1)
            color = func(p0, pe)
            push!(colors, color)

            # face_center = pos .+ [a / 2, (sqrt(3) / 4) * a]
            # @show face_center
            # push!(centers, face_center)

            vi += 1
        end
    end

    faces = []
    for ri in 1:N
        # bulk and right edge
        for j in 2:(N-ri+1-1)
            push!(faces, [vi_rows[ri][j], vi_rows[ri][j+1], vi_rows[ri+1][j]])
            push!(faces, [vi_rows[ri][j], vi_rows[ri+1][j], vi_rows[ri+1][j-1]])
        end
        # left edge
        if ri != N
            push!(faces, [vi_rows[ri][1], vi_rows[ri][2], vi_rows[ri+1][1]])
        end
    end

    length.(vi_rows)
    vi_rows

    vertices = reduce(vcat, transpose.(vertices))
    faces = reduce(vcat, transpose.(faces))

    p = mesh!(ax, vertices, faces;
        color=colors,
        shading=NoShading,
        interpolate=false
    )
    poly!(ax, [(0, 0), (1 / 2, sqrt(3) / 2), (1, 0)];
        color=:transparent,
        strokecolor=:black,
        strokewidth=1.0
    )

    hidedecorations!(ax)

    text!(ax, 0, 0; text="All 0s", align=(:left, :top))
    text!(ax, 1, 0; text="All 1s", align=(:right, :top))
    text!(ax, 1/2, sqrt(3)/2; text="Other acs", align=(:center, :bottom))

    # scatter!(fap.axis, getindex.(centers, 1), getindex.(centers, 2))

    # fap.axis.xticks = range(0., 1., (N-1)+1)
    # fap.axis.yticks = range(0., sqrt(3) / 2, (N-1)+1)

    p
end

function all0sv1s_p1(nmg;
    force_order=true,
    double_arrows=true,
    colors=nothing,
    kwargs...
)
    L = nmg[].L
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    node_clusters = [subgraph(g, (@sprintf "cluster %s" string(i))) for i in 0:L]
    cluster_is = [[] for _ in 0:L]

    for l in labels(nmg)
        label = nmg[l]
        num_ones = count('1', label)
        cluster_index = 1 + num_ones
        push!(cluster_is[cluster_index], label)

        auto_kwargs = (;)
        if !isnothing(colors)
            auto_kwargs = (; auto_kwargs...,
                style="filled",
                color="#" * hex(colors[l]),
                # fontcolor="azure4"
            )
        end

        node_clusters[cluster_index] |> node(label;
            label, auto_kwargs...
        )
    end
    if force_order
        for ci in 1:(length(cluster_is)-1)
            for l1 in cluster_is[ci]
                for l2 in cluster_is[ci+1]
                    g |> edge(l1, l2;
                        style="invis",
                        # weight="100"
                    )
                end
            end
        end
    end

    if double_arrows
        for i in 1:nv(nmg)
            for j in (i+1):nv(nmg)
                il = label_for(nmg, i)
                jl = label_for(nmg, j)
                itoj = haskey(nmg, il, jl)
                jtoi = haskey(nmg, jl, il)
                if itoj && jtoi
                    g |> edge(nmg[il], nmg[jl];
                        dir="both",
                        # weight="10"
                    )
                elseif itoj
                    g |> edge(
                        nmg[il],
                        nmg[jl];
                        # weight="10"
                    )
                elseif jtoi
                    g |> edge(
                        nmg[jl],
                        nmg[il];
                        # weight="10"
                    )
                end
            end
        end
    else
        for (s, d) in edges(graph(ned))
            g |> edge(
                string(e.src),
                string(e.dst);
                # weight="10"
            )
        end
    end

    g
end

function get_all0sv1s_avg_splitprobs(nmg::MetaGraph{C,G,LT};
    cmap0sv1s=ColorSchemes.RdBu_4,
    nondetcolor=ColorSchemes.viridis[end],
) where {C,G,LT}
    L = nmg[].L
    all0sstate = fill(0, L)
    all1sstate = fill(1, L)

    colors = Dict{LT,Color}()
    p0s = [[] for _ in 0:L]
    p1s = [[] for _ in 0:L]
    pxs = [[] for _ in 0:L]

    ff = make_full_splitprober(nmg)

    for l in labels(nmg)
        sps = ff(l)
        prob_all0s = get(sps, [all0sstate], 0.0)
        prob_all1s = get(sps, [all1sstate], 0.0)
        prob_nonextreme = 0.0
        for ac in keys(sps)
            if ac != [all0sstate] && ac != [all1sstate]
                prob_nonextreme += sps[ac]
            end
        end

        num1s = count(x -> x == 1, l)
        push!(p0s[num1s+1], prob_all0s)
        push!(p1s[num1s+1], prob_all1s)
        push!(pxs[num1s+1], prob_nonextreme)
    end

    0:L, mean.(p0s), mean.(p1s), mean.(pxs)
end
