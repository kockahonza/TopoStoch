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
