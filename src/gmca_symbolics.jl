################################################################################
# Extend the functions from gm_symbolics.jl
################################################################################
function get_variables(em::EnergyMatrices)
    union(get_variables(em.monomer), get_variables(em.interactions))
end
function get_variables(ca::ComplexAllosteryGM{S,Num}) where {S}
    variables = Set{Num}()
    for edge in edges(ca.graph)
        for var in Symbolics.get_variables(weight(edge))
            push!(variables, var)
        end
    end
    union(variables, get_variables(ca.energy_matrices))
end

function ssubstitute(em::EnergyMatrices, terms)
    EnergyMatrices(ssubstitute(em.monomer, terms), ssubstitute(em.interactions, terms))
end
function ssubstitute(ca::ComplexAllosteryGM{S,F}, terms::Dict) where {S,F}
    new_energy_matrices = ssubstitute(ca.energy_matrices, terms)
    new_graph = SimpleWeightedDiGraph(
        ssubstitute(adjacency_matrix(ca.graph), terms)
    )
    new_metadata = Dict{Any,Any}(terms)
    new_metadata["old"] = ca.metadata
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(F),
        energy_matrices=new_energy_matrices,
        graph=new_graph,
        metadata=new_metadata
    )
end

function substitute_to_float(em::EnergyMatrices, terms::Dict{Num,Float64})
    EnergyMatrices{Float64}(substitute_to_float(em.monomer, terms), substitute_to_float(em.interactions, terms))
end
function substitute_to_float(ca::ComplexAllosteryGM{S}, terms::Dict{Num,Float64}) where {S}
    new_energy_matrices = substitute_to_float(ca.energy_matrices, terms)
    new_graph = SimpleWeightedDiGraph(
        substitute_to_float(adjacency_matrix(ca.graph), terms)
    )
    new_metadata = Dict{Any,Any}(terms)
    new_metadata["old"] = ca.metadata
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(Float64),
        energy_matrices=new_energy_matrices,
        graph=new_graph,
        environment=(terms[Symbolics.variable(:μ)], terms[Symbolics.variable(:kT)]),
        metadata=new_metadata
    )
end

function make_factory(em::EnergyMatrices, variables=Num.(get_variables(em)); kwargs...)
    mono_fact = make_factory(em.monomer, variables; kwargs...)
    int_fact = make_factory(em.interactions, variables; kwargs...)

    function (args...)
        EnergyMatrices(
            Matrix{Float64}(mono_fact(args...)),
            Matrix{Float64}(int_fact(args...))
        )
    end
end
function make_factory(ca::ComplexAllosteryGM{S}, variables=Num.(get_variables(ca)); kwargs...) where {S}
    em_fact = make_factory(ca.energy_matrices, variables; kwargs...)
    adjmat_fact = make_factory(adjacency_matrix(ca.graph), variables; kwargs...)

    function (args...)
        ComplexAllosteryGM(ca.N, ca.C, ca.B;
            symmetry=S(),
            numtype=Val(Float64),
            energy_matrices=em_fact(args...),
            graph=SimpleWeightedDiGraph(adjmat_fact(args...))
        )
    end
end

################################################################################
# Particular running functions for ComplexAllosteryGM
################################################################################
function ca_eigen1(ca::ComplexAllosteryGM; simplify=false)
    setm, rt = etransmat_safe(ca)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end

function ca_eigen2(ca::ComplexAllosteryGM; simplify=false)
    etm = etransmat(ca)
    setm, rt = substitute_safenames(etm)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end
