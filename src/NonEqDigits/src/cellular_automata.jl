################################################################################
# Basics of relating our Ks to CE numerical codes
################################################################################
function cecode_to_Ks(code)
    if !(0 <= code < 2^8)
        throw(ArgumentError("code must be between 0 and 2^8-1 inclusive"))
    end

    cd = digits(code; base=2, pad=8)

    K01 = [cd[1] cd[2]; cd[5] cd[6]]
    K10 = neg_int.([cd[3] cd[4]; cd[7] cd[8]])

    K01, K10
end
function Ks_to_cecode(K01, K10)
    cd = fill(0, 8)

    cd[1] = K01[1, 1]
    cd[2] = K01[1, 2]
    cd[5] = K01[2, 1]
    cd[6] = K01[2, 2]

    cd[3] = neg_int(K10[1, 1])
    cd[4] = neg_int(K10[1, 2])
    cd[7] = neg_int(K10[2, 1])
    cd[8] = neg_int(K10[2, 2])

    sum(cd .* 2 .^ (0:7))
end
export cecode_to_Ks, Ks_to_cecode

function make_ce_ned(L, code; show=true)
    ned = NonEqDigitsGM(2, L)
    K01, K10 = cecode_to_Ks(code)
    if show
        @show K01
        @show K10
    end
    add_mech_BNN!(ned, Inf, K01)
    add_mech_BNN!(ned, -Inf, K10)
    ned
end
export make_ce_ned

################################################################################
# Simple calculations, stats, properties etc.
################################################################################
function ce_calc_numarrows(K01, K10)
    count(!iszero, K01) + count(!iszero, K10)
end
ce_calc_numarrows(code) = ce_calc_numarrows(cecode_to_Ks(code)...)
export ce_calc_numarrows

################################################################################
# Dealing with symmetries and enumerating all codes
################################################################################
function cesym_01(K01, K10)
    [K10[2, 2] K10[2, 1]; K10[1, 2] K10[1, 1]], [K01[2, 2] K01[2, 1]; K01[1, 2] K01[1, 1]]
end
function cesym_LR(K01, K10)
    transpose(K01), transpose(K10)
end
export cesym_01, cesym_LR

function ce_ned_sym_graph(; noselfsyms=true)
    g = SimpleGraph(2^8)
    for code in 0:2^8-1
        baseKs = cecode_to_Ks(code)
        Ksafter01 = cesym_01(baseKs...)
        KsafterLR = cesym_LR(baseKs...)
        if noselfsyms
            add_edge_ifnotself!(g, code + 1, Ks_to_cecode(Ksafter01...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cecode(KsafterLR...) + 1)
        else
            add_edge!(g, code + 1, Ks_to_cecode(Ksafter01...) + 1)
            add_edge!(g, code + 1, Ks_to_cecode(KsafterLR...) + 1)
        end
    end
    g
end
ce_ned_sym_graphplot(g=ce_ned_sym_graph(); kwargs...) = graphplot(g; nlabels=string.(0:255), kwargs...)
export ce_ned_sym_graph, ce_ned_sym_graphplot

function ce_ucodes_bydeg(by=ce_calc_numarrows)
    comps_by_len = sort(connected_components(ce_ned_sym_graph()); by=length)

    uniques = []

    codes_of_this_len = []
    cur_len = 1
    for cc in comps_by_len
        if length(cc) != cur_len
            sort!(codes_of_this_len; by)
            append!(uniques, codes_of_this_len)
            codes_of_this_len = []
            cur_len = length(cc)
        end
        push!(codes_of_this_len, cc[1] + 1)
    end
    append!(uniques, codes_of_this_len)

    uniques
end
export ce_ucodes_bydeg

################################################################################
# Reading and dealing with the standard CE classifications
################################################################################
const ce_subclasses = ["1", "LP", "HP", "T", "U", "C"]
const ce_class_subclass_pairs = [1 => "1", 2 => "LP", 2 => "HP", 2 => "T", 3 => "C", 3 => "U", 4 => "T", 4 => "C"]
export ce_subclasses, ce_class_subclass_pairs

function ce_load_classes_df()
    DataFrame(CSV.File(datadir("simplececlasses.csv")))
end
export ce_load_classes_df

################################################################################
# Plots
################################################################################
function plot_percecode_vsN(codes, Ns, funcs;
    color=nothing,
    colormap=:managua,
    colorrange=nothing,
    cbarlabel=nothing,
    yoffset_dev=0.0,
    xoffset_dev=0.0,
    do3d=true,
)
    numcodes = length(codes)
    numfuncs = length(funcs)
    # Name unnamed funcs
    namedfuncs = Tuple{String,Function}[]
    func_i = 1
    for func in funcs
        if isa(func, Tuple)
            push!(namedfuncs, func)
        elseif isa(func, Pair)
            push!(namedfuncs, (func[1], func[2]))
        else
            push!(namedfuncs, ((@sprintf "Unnamed %d" func_i), func))
            func_i += 1
        end
    end

    # Prep axes
    (nrows, ncols) = ntogridsize1(numfuncs)
    cis = CartesianIndices((ncols, nrows))
    fig = Figure()
    axs = []
    for ((name, _), ci) in zip(namedfuncs, cis)
        if do3d
            ax = Axis3(fig[ci[2], ci[1]])
        else
            ax = Axis(fig[ci[2], ci[1]])
        end
        ax.title = name
        push!(axs, ax)
    end

    # Prep colormapping
    if isnothing(color)
        color = identity
        if isnothing(colorrange)
            colorrange = (minimum(codes), maximum(codes))
        end
        if isnothing(cbarlabel)
            cbarlabel = "code"
        end
    end
    if isnothing(cbarlabel)
        cbarlabel = "user provided"
    end
    cmap = cgrad(colormap)

    plots = []
    colorvals = []
    for code in codes
        valss = [[] for _ in 1:length(namedfuncs)]
        for N in Ns
            ned = make_ce_ned(N, code; show=false)
            for ((_, func), vals) in zip(namedfuncs, valss)
                push!(vals, func(ned))
            end
        end
        for (ax, vals) in zip(axs, valss)
            cval = color(code)
            if do3d
                ls = lines!(ax, fill(cval, length(Ns)), Ns .+ xoffset_dev * rand(), vals .+ yoffset_dev * rand())
            else
                ls = lines!(ax, Ns .+ xoffset_dev * rand(), vals .+ yoffset_dev * rand())
            end
            push!(plots, ls)
            push!(colorvals, cval)
        end
    end

    if isnothing(colorrange)
        colorrange = (minimum(colorvals), maximum(colorvals))
    end
    colorrange_span = colorrange[2] - colorrange[1]
    colorrange_base = colorrange[1]
    for (p, cv) in zip(plots, colorvals)
        p.color = get(cmap, (cv - colorrange_base) / colorrange_span)
    end

    Colorbar(fig[:, ncols+1]; colorrange, colormap=cmap, label=cbarlabel)

    fig
end
export plot_percecode_vsN
