################################################################################
# Plotting util
################################################################################
struct FigureAxisAnything
    figure::Figure
    axis::Any
    obj::Any
end
display(faa::FigureAxisAnything) = display(faa.figure)

function make_linalpha_cmap(cmap; amin=0.0, amax=1.0)
    cmap = Makie.to_colormap(cmap)
    acmap = [coloralpha(color(c), a) for (c, a) in zip(cmap, LinRange(amin, amax, length(cmap)))]
    cgrad(acmap)
end

function p_get_adjmat_symsafe(graph::AbstractNumGraph)
    map(x -> if iszero(x)
            0.0
        else
            1.0
        end, adjacency_matrix(graph))
end
function p_get_adjmat_symsafe(graph::AbstractFloatGraph)
    adjacency_matrix(graph)
end
p_get_adjmat_symsafe(gm::AbstractGraphModel) = p_get_adjmat_symsafe(graph(gm))

function p_make_ca_ax(dim, place)
    if dim == 2
        Axis(place)
    elseif dim == 3
        LScene(place)
    else
        throw(ErrorException("the layout dimension was not 2 or 3"))
    end
end

function makeprefunclayout(base_layout::NetworkLayout.AbstractLayout, f, args...; kwargs...)
    copyf = copyand(f, args...; kwargs...)
    function (ca)
        base_layout(p_get_adjmat_symsafe(copyf(ca)))
    end
end

# adding extras
function edgecolorbar(fig, plot)
    edgeploti = findfirst(typeof(plot_) <: Plot{GraphMakie.edgeplot} for plot_ in plot.plots)
    Colorbar(fig[1, 2], plot.plots[edgeploti])
end
edgecolorbar((fig, ax, plot)) = edgecolorbar(fig, plot)

function make_interactive(ax, plot)
    interactions_ = keys(interactions(ax))
    if :rectanglezoom in interactions_
        deregister_interaction!(ax, :rectanglezoom)
    end
    if :ndrag in interactions_
        deregister_interaction!(ax, :ndrag)
    end
    function node_drag_action(state, idx, event, axis)
        plot[:node_pos][][idx] = event.data
        plot[:node_pos][] = plot[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)
end
make_interactive((_, ax, plot)) = make_interactive(ax, plot)

################################################################################
# Main plotting
################################################################################

################################################################################
# Util
################################################################################
function savefig(subdir, prefix, gm::AbstractGraphModel, fig)
    savefig(subdir,
        savename(prefix, gm;
            allowedtypes=[DrWatson.default_allowed(gm)...; Dict],
            expand=["metadata"]),
        fig)
end
