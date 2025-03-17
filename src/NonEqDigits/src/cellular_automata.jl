################################################################################
# Basics of relating our Ks to CA numerical codes
################################################################################
function cacode_to_Ks(code)
    if !(0 <= code < 2^8)
        throw(ArgumentError("code must be between 0 and 2^8-1 inclusive"))
    end

    cd = digits(code; base=2, pad=8)

    K01 = [cd[1] cd[2]; cd[5] cd[6]]
    K10 = neg_int.([cd[3] cd[4]; cd[7] cd[8]])

    K01, K10
end
function Ks_to_cacode(K01, K10)
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
export cacode_to_Ks, Ks_to_cacode

function make_ca_ned(L, code; show=true)
    ned = NonEqDigitsGM(2, L)
    K01, K10 = cacode_to_Ks(code)
    if show
        @show K01
        @show K10
    end
    add_mech_BNN!(ned, Inf, K01)
    add_mech_BNN!(ned, -Inf, K10)
    ned
end
export make_ca_ned

################################################################################
# Simple calculations, stats, properties etc.
################################################################################
function ca_calc_numarrows(K01, K10)
    count(!iszero, K01) + count(!iszero, K10)
end
ca_calc_numarrows(code) = ca_calc_numarrows(cacode_to_Ks(code)...)
export ca_calc_numarrows

################################################################################
# Dealing with symmetries and enumerating all codes
################################################################################
function caKssym_01(K01, K10)
    [K10[2, 2] K10[2, 1]; K10[1, 2] K10[1, 1]], [K01[2, 2] K01[2, 1]; K01[1, 2] K01[1, 1]]
end
function caKssym_LR(K01, K10)
    transpose(K01), transpose(K10)
end
function caKssym_F(K01, K10)
    caKssym_LR(caKssym_01(K01, K10)...)
end
export caKssym_01, caKssym_LR, caKssym_F

# the extra, time symmetries coming into ca_ned_ultrasym_graph
caKssym_T(K01, K10) = K10, K01
caKssym_T01(args...) = caKssym_T(caKssym_01(args...)...)
caKssym_TLR(args...) = caKssym_T(caKssym_LR(args...)...)
caKssym_TF(args...) = caKssym_T(caKssym_F(args...)...)
export caKssym_T, caKssym_T01, caKssym_TLR, caKssym_TF

function ca_has01sym(args...)
    args == caKssym_01(args...)
end
ca_has01sym(code) = ca_has01sym(cacode_to_Ks(code)...)
function ca_hasLRsym(args...)
    args == caKssym_LR(args...)
end
ca_hasLRsym(code) = ca_hasLRsym(cacode_to_Ks(code)...)
function ca_hasFsym(args...)
    args == caKssym_F(args...)
end
ca_hasFsym(code) = ca_hasFsym(cacode_to_Ks(code)...)
ca_hasallsyms(args...) = ca_hasLRsym(args...) && ca_has01sym(args...)
export ca_has01sym, ca_hasLRsym, ca_hasFsym, ca_hasallsyms

function ca_symclass(args...)
    has01sym = ca_has01sym(args...)
    hasLRsym = ca_hasLRsym(args...)
    if has01sym && hasLRsym
        :all
    elseif hasLRsym
        :h
    elseif has01sym
        :g
    elseif ca_hasFsym(args...)
        :f
    else
        :none
    end
end
export ca_symclass

ca_iseq(K01, K10) = K01 == K10
ca_iseq(code) = ca_iseq(cacode_to_Ks(code)...)
export ca_iseq

ca_numenzymes(Ks::Vararg{Any,2}) = sum(count.(x -> x != 0, Ks))
ca_numenzymes(code) = ca_numenzymes(cacode_to_Ks(code)...)
export ca_numenzymes

function ca_ned_sym_graph(; noselfsyms=true)
    g = SimpleGraph(2^8)
    for code in 0:2^8-1
        baseKs = cacode_to_Ks(code)
        Ksafter01 = caKssym_01(baseKs...)
        KsafterLR = caKssym_LR(baseKs...)
        KsafterF = caKssym_F(baseKs...)
        if noselfsyms
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(Ksafter01...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterLR...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterF...) + 1)
        else
            add_edge!(g, code + 1, Ks_to_cacode(Ksafter01...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterLR...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterF...) + 1)
        end
    end
    g
end
ca_ned_sym_graphplot(g=ca_ned_sym_graph(); kwargs...) = graphplot(g; nlabels=string.(0:255), kwargs...)
export ca_ned_sym_graph, ca_ned_sym_graphplot

function ca_ucodes_sorted(by=identity)
    sort(minimum.(connected_components(ca_ned_sym_graph())) .- 1; by)
end
function ca_ucodes_bydeg(by=ca_calc_numarrows)
    comps_by_len = sort(connected_components(ca_ned_sym_graph()); by=length)

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
        # always take the smallest code from each symmetry component
        push!(codes_of_this_len, minimum(cc) - 1)
    end
    append!(uniques, codes_of_this_len)

    uniques
end
export ca_ucodes_sorted, ca_ucodes_bydeg

function ca_ned_ultrasym_graph(; noselfsyms=true)
    g = SimpleGraph(2^8)
    for code in 0:2^8-1
        baseKs = cacode_to_Ks(code)
        Ksafter01 = caKssym_01(baseKs...)
        KsafterLR = caKssym_LR(baseKs...)
        KsafterF = caKssym_F(baseKs...)
        KsafterT01 = caKssym_T01(baseKs...)
        KsafterTLR = caKssym_TLR(baseKs...)
        KsafterTF = caKssym_TF(baseKs...)
        if noselfsyms
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(Ksafter01...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterLR...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterF...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterT01...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterTLR...) + 1)
            add_edge_ifnotself!(g, code + 1, Ks_to_cacode(KsafterTF...) + 1)
        else
            add_edge!(g, code + 1, Ks_to_cacode(Ksafter01...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterLR...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterF...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterT01...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterTLR...) + 1)
            add_edge!(g, code + 1, Ks_to_cacode(KsafterTF...) + 1)
        end
    end
    g
end
ca_ned_ultrasym_graphplot(g=ca_ned_ultrasym_graph(); kwargs...) = graphplot(g; nlabels=string.(0:255), kwargs...)
export ca_ned_ultrasym_graph, ca_ned_ultrasym_graphplot

function ca_ultraucodes_sorted(by=identity)
    sort(minimum.(connected_components(ca_ned_sym_graph())) .- 1; by)
end
export ca_ultraucodes_sorted

module CANEDSymGroup
using ..NonEqDigits

function sm01()
    mat = fill(0, 8, 8)
    for i in 1:8
        mat[i, 9-i] = 1
    end
    mat
end
function smLR()
    mat = fill(0, 8, 8)
    mat[1, 1] = 1
    mat[2, 3] = 1
    mat[3, 2] = 1
    mat[4, 4] = 1
    mat[5, 5] = 1
    mat[6, 7] = 1
    mat[7, 6] = 1
    mat[8, 8] = 1
    mat
end
export sm01, smLR

smG = sm01
smH = smLR
function smF()
    sm01() * smLR()
end
export smG, smH, smF

function smT()
    # hcat(vcat(fill(0,4,4),I(4)), vcat(I(4), fill(0,4,4)))
    mat = fill(0, 8, 8)
    for i in 1:4
        mat[i+4, i] = 1
        mat[i, i+4] = 1
    end
    mat
end
export smT

# Can confirm by making a multiplication table that this is in fact what
# I thought it was, aka a finite abelian group with 3 non-identity elements
# each with rank 2. (didn't check associativity tbf)

end

################################################################################
# Reading and dealing with the standard CA classifications
################################################################################
const ca_subclasses = ["1", "LP", "HP", "T", "U", "C"]
const ca_class_subclass_pairs = [1 => "1", 2 => "LP", 2 => "HP", 2 => "T", 3 => "C", 3 => "U", 4 => "T", 4 => "C"]
export ca_subclasses, ca_class_subclass_pairs

function ca_load_classes_df()
    DataFrame(CSV.File(datadir("simplecaclasses.csv")))
end
export ca_load_classes_df

################################################################################
# Plots
################################################################################
function plot_eachcode_vsN(codes, Ns, funcs;
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
            ned = make_ca_ned(N, code; show=false)
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
export plot_eachcode_vsN

function plot_eachcode_vsN2(codes, Ns, funcs;
    yoffset_dev=0.0,
    xoffset_dev=0.0,
)
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
    (nrows, ncols) = ntogridsize1(length(funcs))
    cis = CartesianIndices((ncols, nrows))
    fig = Figure()
    axs = []
    for ((name, _), ci) in zip(namedfuncs, cis)
        ax = Axis(fig[ci[2], ci[1]])
        ax.title = name
        push!(axs, ax)
    end

    plots = []
    for (code_i, code) in enumerate(codes)
        valss = [[] for _ in 1:length(namedfuncs)]
        for N in Ns
            ned = make_ca_ned(N, code; show=false)
            for ((_, func), vals) in zip(namedfuncs, valss)
                push!(vals, func(ned))
            end
        end
        for (ax, vals) in zip(axs, valss)
            ls = lines!(ax, Ns .+ xoffset_dev * (code_i - 1), vals .+ yoffset_dev * (code_i - 1); label=f"{code}")
            push!(plots, ls)
        end
    end

    Legend(fig[:, ncols+1], axs[1], "Rules")

    fig
end
export plot_eachcode_vsN2
