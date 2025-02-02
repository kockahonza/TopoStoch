using DrWatson
@quickactivate "TopoStochSim"

using CairoMakie
using GLMakie
using ColorBrewer # For single colors

using ComplexAllostery
using ComplexAllostery.Model2_5_C2

################################################################################
# 24-12-03 - looking at the cycles in N~3-6 and the pertrubation states
function m1()
    ca = makeCAGM(4, 1; vars=vars_simplified(; newrs=:ss1, energies=false, concentrations=:ss1))
    @show get_variables(ca)
    cca = substitute_to_float(ca, [0, 0, 2, 10, 0.001])

    sc1 = GLMakie.Screen()
    display(sc1, plotcgm(cca; layout=:NTc))

    es = fulleigenanalysis(cca)[2]
    @show findall(isnothing, degenerate_map(es))

    display(plot_ps(cca, evec(es, 1); layout=:NTc))

    ca, cca, es, sc1
end
# ca, cca, es, sc1 = m1()

function ap1(N=4; kwargs...)
    ca = makeCAGM(N, 1; vars=vars_simplified(; newrs=:ss1, energies=false, concentrations=:ss1))

    @show get_variables(ca)
    cca = substitute_to_float(ca, [0, 0, 3, 10, 0.0])
    ss = supersteadystate(cca)

    # Plotting starts
    fig = Figure()
    diagram_ax = Axis(fig[2, 1]; aspect=1.0)
    data_ax = p_make_ax(2, fig[:, 2])
    hidedecorations!(diagram_ax)
    hidedecorations!(data_ax)

    colsize!(fig.layout, 1, Relative(0.5))
    rowsize!(fig.layout, 2, Aspect(1, 1.0))

    ############################################################
    # Make the diagram
    fakeN = 6
    monomer_size = 1.0
    polymer_rad = 1.1
    ligand_size = 0.2

    pal = palette("Dark2", 5)
    ligand_color = pal[3]
    relaxed_color = pal[1]

    base_attrs = (;
        color=nothing,
        strokecolor=:black,
        strokewidth=1,
    )
    function Tat(center; kwargs...)
        poly!(diagram_ax,
            Rect2(center .- (monomer_size / 2), monomer_size, monomer_size);
            base_attrs...,
            kwargs...
        )
    end
    function Rat(center; kwargs...)
        poly!(diagram_ax,
            Circle(center, 1.05 * monomer_size / 2);
            base_attrs...,
            color=alphacolor(relaxed_color, 0.3),
            kwargs...
        )
    end
    function Lat(center; kwargs...)
        poly!(diagram_ax,
            Circle(center, ligand_size / 2);
            base_attrs...,
            color=alphacolor(ligand_color, 0.5),
            kwargs...
        )
    end

    # calc positions
    positions = Vector{Point{2,Float32}}(undef, fakeN)
    ligand_positions = Vector{Point{2,Float32}}(undef, fakeN)
    for i in 1:fakeN
        th = ((i - 1) / fakeN) * 2 * pi
        positions[i] = polymer_rad .* Point2f(cos(th), sin(th))
        ligand_positions[i] = positions[i] + ((monomer_size / 2) + 0.4 * ligand_size) .* Point2f(cos(th), sin(th))
    end

    # draw monomers
    positionsT = positions[begin+3:end]
    positionsR = positions[begin:begin+2]
    map(Tat, positionsT)
    map(Rat, positionsR)

    # draw faint line connecting monomers
    for i in 1:fakeN
        lines!(diagram_ax, [positions[i], positions[mod1(i + 1, fakeN)]],
            color=:gray,
            linestyle=:dash
        )
    end

    # draw ligands
    positionsL = ligand_positions[[1, 3, 4]]
    map(Lat, positionsL)

    # draw labels
    diagram_ax.title = "Example assembly microstate\n(corresponds to a 3,3 state)"
    position_Rlabel = Point2f(1.2, 2.0)
    position_Tlabel = Point2f(-1.2, -2.0)

    text!(diagram_ax, position_Rlabel; text="R conformation", align=(:center, :bottom))
    dirs = [posR .- position_Rlabel for posR in positionsR]
    arrows!(diagram_ax, fill(position_Rlabel, length(positionsR)), 0.8 .* dirs; color=:gray,)

    text!(diagram_ax, position_Tlabel; text="T conformation", align=(:center, :top))
    dirs = [posT .- position_Tlabel for posT in positionsT]
    arrows!(diagram_ax, fill(position_Tlabel, length(positionsT)), 0.8 .* dirs; color=:gray,)

    position_Llabel = Point2f(-1.1, 2.0)
    text!(diagram_ax, position_Llabel; text="bound ligands", align=(:center, :bottom))
    dirs = [posL .- position_Llabel for posL in positionsL]
    arrows!(diagram_ax, fill(position_Llabel, length(positionsL)), dirs .- ligand_size * normalize.(dirs); color=:gray,)


    # final bits
    xlims!(diagram_ax, -2.5, 2.5)
    ylims!(diagram_ax, -2.5, 2.5)

    ############################################################
    # Make the data plot
    function remove_bulk_labels(akwargs)
        labels = akwargs[:nlabels]
        for (i, st) in enumerate(allstates(cca))
            numlig = calc_numligands(st)
            numR = calc_numofconf(st, 2)
            if (numlig != 0) && (numlig != N) && (numR != 0) && (numR != N)
                labels[i] = ""
            end
        end
    end
    cplot = plotcgm!(data_ax, cca, ss;
        layout=Spectral(dim=2),
        node_size=10.0,
        nlabels_align=(:center, :center),
        nlabels_color=[[ligand_color, :black, :black, relaxed_color] for _ in 1:numstates(cca)],
        makwargs=remove_bulk_labels,
        kwargs...
    )

    # do reasonable offsets
    offsets = [Point2f(0.0, 0.0) for _ in 1:numstates(cca)]
    labeled_already = []
    for (i, st) in enumerate(allstates(cca))
        numlig = calc_numligands(st)
        numR = calc_numofconf(st, 2)
        if !((numlig != 0) && (numlig != N) && (numR != 0) && (numR != N))
            macrostate = (numlig, numR)
            if !(macrostate in labeled_already)
                node_pos = cplot.node_pos[][i]
                offsets[i] = 0.08 * normalize(node_pos)
                push!(labeled_already, macrostate)
            else
                cplot.nlabels[][i] = ""
            end
        end
    end
    cplot.nlabels_offset = offsets
    cplot.nlabels = cplot.nlabels[]

    data_ax.title = f"Steady state probability current network between\nall (overlaid) microstates with {N} subunits"

    ############################################################
    # final bits
    FigureAxisAnything(fig, diagram_ax, [diagram_ax, data_ax, cplot])
end
