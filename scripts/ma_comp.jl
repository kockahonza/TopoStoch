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
        graph_data="ma comp graph"
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
        reduce_allcodes_compgraph(cgraph)
    end
end

function reduce_allcodes_compgraph(accg)
    scgraph = MetaGraph(
        DiGraph();
        label_type=Int,
        vertex_data_type=String,
        edge_data_type=Tuple{Int,Float64},
        graph_data="ma comp graph"
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
function condense_acs(g)
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
