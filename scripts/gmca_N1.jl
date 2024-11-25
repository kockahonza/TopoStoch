using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(scriptsdir("gmca_v2_5.jl"))

function N1sscur(; dofull=false)
    ca, sca = make_v2_5(1, 1; simplified=:both)
    tm, rt = etransmat_safe(ca; mat=true)
    tm2, rt2 = etransmat_safe(sca; mat=true)
    @assert isequal(tm, tm2)

    @variables a, b, c, d, e, f, g, h
    es = w_eigen(tm; safe=false)
    @assert es.evals[1] == 0
    ssv = es.evecs[1] * (a * d * e + a * d * f + b * c * f + b * d * f)
    ssn = w_simplify(sum(ssv); safe=false)

    cs = simplify.(calc_currents(tm, ssv); expand=true)
    sscv = w_simplify(ssubstitute(cs[3, 1], rt))
    if dofull
        sscn = w_simplify(ssubstitute(ssn, rt))
    else
        sscn = nothing
    end

    sscv2 = w_simplify(ssubstitute(cs[3, 1], rt2))
    sscn2 = w_simplify(ssubstitute(ssn, rt2))

    (; ca, sca, tm, rt, tm2, rt2, es, ssv, ssn, cs, sscv, sscn, sscv2, sscn2)
end

function N1eqss((; ssv, ssn, rt2))
    essv = w_simplify(ssubstitute(w_simplify(ssubstitute(ssv, rt2); full=true), get_eq_subs1()); full=true)
    essn = w_simplify(ssubstitute(w_simplify(ssubstitute(ssn, rt2); full=true), get_eq_subs1()); full=true)

    ess = w_simplify(essv / essn; full=true)

    (; essv, essn, ess)
end
N1eqss(_::Nothing) = eqss(N1sscur())

function N1curplot(type=:lines;
    crange=LinRange(0.0, 10.0, 100),
    keeprs=true,
    colormap=:blues,
    returnsgo=false
)
    aa = N1sscur()

    sscur = aa.sscv2 / aa.sscn2
    sscur = ssubstitute(sscur, get_concentrations()[3] => 0.0)
    vars = get_sem_vars()[1:2]
    ranges = [0.0:0.1:10.0, 0.0:0.05:5.0]
    startvalues = [0.0, 0.0]
    if !keeprs
        sscur = ssubstitute(sscur, ssr_terms())
    else
        vars = [get_newrs(); vars]
        ranges = [[0.0:0.01:3.0, 0.0:0.01:3.0, 0.0:0.01:3.0, 0.0:0.01:2.0]; ranges]
        startvalues = [[1.0, 1.0, 1.0, 0.0]; startvalues]
    end
    sscur = w_simplify(sscur)

    conc_vars = get_concentrations()[1:2]
    sscurf = make_factory(sscur, [vars; conc_vars])

    fig, sgo, sg = prep_sliders(vars; ranges, startvalues)

    if type in SA[:contourf, :heatmap, :surface]
        sscur_vals = lift(sgo...) do (svals...)
            sscurf.(svals..., crange, transpose(crange))
        end
        curmax = maximum(sscur_vals[])
    end

    if type == :contourf
        ax = p_make_ax(2, fig[2, 1], string.(reverse(conc_vars)))
        obj = contourf!(ax, crange, crange, lift(transpose, sscur_vals);
            levels=(LinRange(-curmax, curmax, 40)),
            colormap=:redsblues
        )

        Colorbar(fig[2, 2], obj)
    elseif type == :heatmap
        ax = p_make_ax(2, fig[2, 1], string.(reverse(conc_vars)))
        obj = heatmap!(ax, crange, crange, lift(transpose, sscur_vals);
            colorrange=(-curmax, curmax),
            colormap=:redsblues
        )

        Colorbar(fig[2, 2], obj)
    elseif type == :surface
        ax = p_make_ax(3, fig[2, 1], string.([reverse(conc_vars); "j"]); useaxis3=true)
        obj = surface!(ax, crange, crange, lift(transpose, sscur_vals);
            colorrange=(-curmax, curmax),
            colormap=:redsblues
        )
        cspan = crange[end] - crange[begin]
        # scale!(ax.scene, 1 / cspan, 1 / cspan, 2 * curmax)

        Colorbar(fig[2, 2], obj)
    elseif type == :lines
        ax = p_make_ax(2, fig[2, 1], ["cATP", "current"])
        obj = []

        cmap = PlotUtils.get_colorscheme(colormap)
        cmin = minimum(crange)
        cmax = maximum(crange)

        for cP in crange
            sscur_at_cP = lift(sgo...) do (svals...)
                sscurf.(svals..., cP, crange)
            end
            color = get(cmap, (cP - cmin) / (cmax - cmin))
            push!(obj, lines!(ax, crange, sscur_at_cP; color))
        end

        sscur_at_1 = lift(sgo...) do (svals...)
            sscurf.(svals..., 1.0, crange)
        end
        push!(obj, lines!(ax, crange, sscur_at_1; color=:red))

        Colorbar(fig[2, 2]; label="cP", colormap=cmap, colorrange=(cmin, cmax))
        text!(ax, cmax * 0.95, lift(x -> x[end], sscur_at_1); text="cP=1", color=:red)
    end

    if !returnsgo
        FigureAxisAnything(fig, ax, obj)
    else
        FigureAxisAnything(fig, ax, obj), sgo
    end
end

function make_sscurfactory()
    aa = N1sscur()
    sscur = aa.sscv2 / aa.sscn2

    # substitute cADP = 0
    sscur = w_simplify(ssubstitute(sscur, get_concentrations()[3] => 0.0))

    # vars: r1, r2, r3, alpha, et, der
    vars = [get_newrs(); get_sem_vars()[1:2]]

    # conc_vars eP, eATP
    conc_vars = get_concentrations()[1:2]
    sscurf = make_factory(sscur, [vars; conc_vars])
end

function makeN1curplots1(;
    crange=LinRange(0.0, 10.0, 100),
    et=0.0,
    der=0.0,
    colormap=:blues
)
    sscur_fac = make_sscurfactory()

    params = []
    for alpha in [0.0, 0.1, 0.4]
        for r1 in [1.0, 2.0]
            for r2 in [1.0]
                for r3 in [1.0, 2.0]
                    push!(params, [r1, r2, r3, alpha, et, der])
                end
            end
        end
    end

    nplots = length(params)
    nrows = round(Int, sqrt(nplots), RoundUp)

    cmap = PlotUtils.get_colorscheme(colormap)
    cmin = minimum(crange)
    cmax = maximum(crange)

    fig = Figure()
    for (i, ps) in enumerate(params)
        (r, c) = fldmod1(i, nrows)

        ax = p_make_ax(2, fig[r, c], ["cATP", "current"])
        ax.title = join(ps[begin:end-2], ", ")

        for cP in crange
            color = get(cmap, (cP - cmin) / (cmax - cmin))
            lines!(ax, crange, sscur_fac.(ps..., cP, crange); color)
        end

        sscurr_at_1 = sscur_fac.(ps..., 1.0, crange)
        lines!(ax, crange, sscurr_at_1; color=:red)

        # text!(ax, cmax * 0.95, sscurr_at_1[end]; text="cP=1", color=:red)
    end
    Colorbar(fig[:, div(nplots, nrows)+2]; label="cP", colormap=cmap, colorrange=(cmin, cmax))

    fig
end
