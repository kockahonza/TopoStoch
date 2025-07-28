using DrWatson;
@quickactivate "TopoStochSim";

using GraphModels
using NonEqDigits
using GraphvizDotLang: digraph, subgraph, node, edge, attr
using ColorSchemes

function plot_0s_v_1s_gv(ned::NonEqDigitsGM{Loop,2,L};
    force_order=true,
    double_arrows=true,
    kwargs...
) where {L}
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    node_clusters = [subgraph(g, (@sprintf "cluster %s" string(i))) for i in 0:L]
    cluster_is = [[] for _ in 0:L]

    for i in 1:numstates(ned)
        st = itostate(i, ned)
        num_ones = count(x -> x == 1, st)
        cluster_index = 1 + num_ones
        node_clusters[cluster_index] |> node(string(i); label=join(st))
        push!(cluster_is[cluster_index], i)
    end
    if force_order
        for ci in 1:(length(cluster_is)-1)
            for i1 in cluster_is[ci]
                for i2 in cluster_is[ci+1]
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
        for e in edges(graph(ned))
            g |> edge(
                string(e.src),
                string(e.dst);
                # weight="10"
            )
        end
    end

    g
end

function molaut_gv(ma::MolAut;
    cluster=false,
    force_layout=false,
    simple_nodes=false,
    highlight_acs=nothing,
    edge_colormap=nothing,
    kwargs...
)
    mg = ma.mg

    g = digraph(;
        rankdir="LR",
        # ranksep="1",
        kwargs...
    )

    if simple_nodes
        g |> attr(:node;
            shape="point",
            fixedsize="true",
            width="0.1",
        )
    end

    # prep node styling/coloring etc
    if highlight_acs == true
        highlight_acs = "crimson"
    end
    if !isnothing(highlight_acs)
        acs = attracting_components(mg)
        vs_in_acs = reduce(vcat, acs)
    end
    function node_style(l)
        rslt = Dict{Symbol,Any}()

        if !isnothing(highlight_acs)
            if code_for(mg, l) in vs_in_acs
                rslt[:color] = highlight_acs
            end
        end

        rslt
    end

    # add all nodes to graph
    if cluster
        clusters = [subgraph(g, "cluster_$(i)_ones") for i in 0:ma.L]

        if force_layout
            fake_cluster_nodes = [(@sprintf "fcnode %d" i) for i in 0:ma.L]
            for (fcn, c) in zip(fake_cluster_nodes, clusters)
                c |> node(fcn;
                    shape="point",
                    style="invis"
                )
            end
            for i in 2:length(fake_cluster_nodes)
                g |> edge(fake_cluster_nodes[i-1], fake_cluster_nodes[i];
                    style="invis",
                    weight="100"
                )
            end
        end

        for l in labels(mg)
            num1s = count(x -> x == 1, l)
            c = clusters[1+num1s]
            c |> node(mg[l];
                node_style(l)...
            )
        end
    else
        for l in labels(mg)
            g |> node(mg[l];
                node_style(l)...
            )
        end
    end

    if !isnothing(edge_colormap)
        if edge_colormap isa Symbol
            edge_colormap = colorschemes[edge_colormap]
        end
        edge_vals = unique(Graphs.weights(mg))
        mn, mx = extrema(edge_vals)
        if isapprox(mn, mx)
            @warn "skipping coloring as all edged have the same value"
            cfunc = nothing
        else
            cfunc = cval -> "#" * hex(get(edge_colormap, (cval - mn) / (mx - mn)))
        end
    else
        cfunc = nothing
    end

    # prep edge styling
    function edge_style(s, d)
        rslt = Dict{Symbol,Any}()

        if !isnothing(cfunc)
            rslt[:color] = cfunc(mg[s, d])
        end

        if !isnothing(highlight_acs)
            if (code_for(mg, s) in vs_in_acs) && (code_for(mg, d) in vs_in_acs)
                rslt[:color] = highlight_acs
            end
        end

        rslt
    end

    # add edges to graph
    for (s, d) in edge_labels(mg)
        g |> edge(mg[s], mg[d];
            edge_style(s, d)...
        )
    end

    g
end
