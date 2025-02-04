using DrWatson
@quickactivate "TopoStochSim"

using JLD2
using GLMakie
using Printf

using GraphModels # not needed, just here to help the LSP
using ComplexAllostery
using ComplexAllostery.Model3
using LetterCodes

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

    fap = heatmap(tm; colormap=Makie.Categorical(:tab20), colorrange=(1, 20), axis=(; aspect=1))
    Colorbar(fap.figure[1, 2], fap.plot)
    fap
    # scatter(1:10, 1:10; color=1:10, colormap=:tab20, colorrange=(1,20))
end
