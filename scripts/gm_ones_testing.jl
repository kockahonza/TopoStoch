using DrWatson
@quickactivate "TopoStochSim"

using Ones

# Making OneGMs
function make_big(N)
    ogm = OnesGM(N; metadata=Dict("ctype" => "full"))
    add_edges_cycle!(ogm, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    add_edges_cycle!(ogm, [[1, 0, 0], [1, 1, 0], [1, 0, 1]])
    add_edges_cycle!(ogm, [[0, 1, 0], [1, 1, 0], [0, 1, 1]])
    add_edges_cycle!(ogm, [[0, 0, 1], [1, 0, 1], [0, 1, 1]])
    add_edges_cycle!(ogm, [[0, 1, 1], [1, 1, 1]])
    add_edges_cycle!(ogm, [[1, 0, 1], [1, 1, 1]])
    add_edges_cycle!(ogm, [[1, 1, 0], [1, 1, 1]])
    ogm
end

function make_3denzym(N, args...; kwargs...)
    cycle = [
        [1, 0, 0],
        [1, 1, 0],
        [1, 1, 1],
        [1, 0, 1]
    ]
    make_single_cycle(N, cycle, args...; kwargs...)
end

function make_3denzym2(N, symmetry, alpha; kwargs...)
    cycle1 = [
        [1, 0, 0],
        [1, 1, 0],
        [1, 1, 1],
        [1, 0, 1]
    ]
    cycle2 = [
        [0, 0, 0],
        [0, 1, 0],
        [0, 1, 1],
        [0, 0, 1]
    ]
    make_multi_cycle(N, symmetry, cycle1, 1.0, cycle2, alpha)
end

# Plotting
function makespringplots3()
    for i in 1:4
        for n in 3:6
            ogm = make_single_cycle(n, get_c2d_int(i))
            fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr, layout=:tree, n_ss_colorbar=false)
            savefig(f"ones/tc{i}/", "spring", ogm, fap.figure)
        end
    end
end
