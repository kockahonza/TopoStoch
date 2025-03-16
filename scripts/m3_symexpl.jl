using DrWatson
@quickactivate "TopoStochSim"

using JLD2
using Printf

using GraphModels # not needed, just here to help the LSP
using ComplexAllostery
using ComplexAllostery.Model3
using LetterCodes

if !(@isdefined nomakie) || !nomakie
    println("A")
    using GLMakie
end

################################################################################
# 25-02-04 - looking at the symbolic transmat for patterns
function getca()
    ca = makeCAGM(4, vars_simplified())
    @show get_variables(ca)
    ca
end

function unique_sym_transmat(letters=false; mat=true)
    ca = getca()

    tm = transmat(ca; mat)
    uniques = unique(tm)
    @printf "num of unique is %d, should be 8\n" length(uniques)

    terms = Dict{Num,Num}()
    if letters
        code = Code()
        for uniq in uniques
            if !iszero(uniq)
                terms[uniq] = Symbolics.variable(string(code))
            end
            increment!(code)
        end
    else
        num = 1
        for uniq in uniques
            if !iszero(uniq)
                terms[uniq] = num
            end
            num += 1
        end
    end

    ssubstitute(tm, terms)
end

function uniqtmheatmap(halfplane=true)
    ca = getca()
    tm = transmat(ca; mat=true)

    uniques = unique(tm)
    # terms = Dict{Num,Num}(uniques .=> [NaN, -3, +1, +3, 4, +2, -1, -2])
    terms = Dict{Num,Num}(uniques .=> [NaN, 6, 1, 5, 7, 3, 2, 4])
    display(terms)

    tm = ssubstitute(tm, terms)
    for x in 1:size(tm)[1]
        for y in 1:size(tm)[2]
            if y > x
                tm[x, y] = NaN
            end
        end
    end

    colormap = cgrad(colors_tab20[begin:begin+6], 7; categorical=true)
    fap = heatmap(tm; colormap, axis=(; aspect=1))
    cb = Colorbar(fap.figure[1, 2], fap.plot)

    # adjust the colorbar to be fake categorical
    cb.axis.attributes.limits = (0.5, 7.5)
    color_labels = [
        "1 - +P,tense",
        "2 - -P,tense",
        "3 - +P,relaxed",
        "4 - -P,relaxed",
        "5 - C,out of bulk",
        "6 - C,into bulk",
        "7 - C,within bulk",
    ]
    println(color_labels)
    cb.ticks = (1:7, color_labels)

    fap
    # To test the colors are working correctly
    # scatter(1:10, 1:10; color=1:10, colormap=:tab20, colorrange=(1,20))
end

################################################################################
# 25-02-05 - just basic exploration
function getcaandfact(N=4)
    ca = makeCAGM(N, vars_simplified())
    @show get_variables(ca)

    ordered_vars = [symvars.sim_rs; symvars.energies[3]; symvars.concentrations[1]; symvars.redconcentration]
    @show ordered_vars

    ca, make_factory(ca, ordered_vars)
end
