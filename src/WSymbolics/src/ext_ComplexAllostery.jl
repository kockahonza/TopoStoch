using .ComplexAllostery

function ca_eigen1(ca::ComplexAllosteryGM; simplify=false)
    setm, rt = etransmat_safe(ca)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end
export ca_eigen1

function ca_eigen2(ca::ComplexAllosteryGM; simplify=false)
    etm = etransmat(ca)
    setm, rt = substitute_safenames(etm)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end
export ca_eigen2
