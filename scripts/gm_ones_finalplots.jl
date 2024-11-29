using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gm_ones.jl"))

function make_c2d_simplecyle_plots()
    for i in 1:4
        for n in 3:8
            ogm = make_single_cycle(n, get_c2d_int(i))
            fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr, n_ss_colorbar=false)
            savefig(f"ones/c{i}/", "spring", ogm, fap.figure)
        end
    end
end

function make_c2d_simplecyle_chain_plots()
    for i in 1:4
        for n in 3:8
            ogm = make_single_cycle(n, get_c2d_int(i), Chain())
            fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr, n_ss_colorbar=false)
            savefig(f"ones/cc{i}/", "spring", ogm, fap.figure)
        end
    end
end

function make_c2d_simplecyle_weighted_plots()
    for n in 3:6
        for i in 1:4
            wcycle = [(c, 1.0 + 0.5 * Int(ci == i)) for (ci, c) in enumerate(get_c2d_int(1))]
            ogm = make_single_wcycle(n, wcycle)
            fap = plotgm(ogm; plotgm_save_kwargs()..., e_color=:rates, layout=:tree)
            savefig("ones/c1weights/", f"R_N={ogm.N}w={i}", fap.figure)
            fap = plotgm(ogm; plotgm_save_kwargs()..., e_color=:currents, layout=:tree)
            savefig("ones/c1weights/", f"N={ogm.N}w={i}", fap.figure)
        end
    end
end

function main()
    make_c2d_simplecyle_plots()
    make_c2d_simplecyle_chain_plots()
    make_c2d_simplecyle_weighted_plots()
end
