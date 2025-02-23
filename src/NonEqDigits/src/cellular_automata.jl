function cecode_to_2Ks(code)
    if !(0 <= code < 2^8)
        throw(ArgumentError("code must be between 0 and 2^8-1 inclusive"))
    end

    cd = digits(code; base=2, pad=8)

    K01 = [cd[1] cd[2]; cd[5] cd[6]]
    K10 = neg_int.([cd[3] cd[4]; cd[7] cd[8]])

    K01, K10
end
export cecode_to_2Ks

function make_ce_ned(L, code)
    ned = NonEqDigitsGM(2, L)
    K01, K10 = cecode_to_2Ks(code)
    @show K01
    @show K10
    add_mech_BNN!(ned, Inf, K01)
    add_mech_BNN!(ned, -Inf, K10)
    ned
end
export make_ce_ned
