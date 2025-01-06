
mutable struct CellSystem{D, LU, TB, TC} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    cell::TC
    function CellSystem(sys::AbstractIsolatedSystem{D, LU}, cell::PeriodicCell{D,T}) where {D, LU, T}
        new{D, LU, typeof(sys), typeof(cell)}(sys, cell)
    end
end

function CellSystem(sys::Union{AbstractSystem,AtomsVector}, i)
    base_sys = AtomicPropertySystem(sys, i)
    if isa(cell(sys), IsolatedCell)
        return base_sys
    end
    return CellSystem(base_sys, cell(sys))
end

function CellSystem(sys::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies)
    base_sys = AtomicPropertySystem(sys, spc)
    if isa(cell(sys), IsolatedCell)
        return base_sys
    end
    return CellSystem(base_sys, cell(sys))
end

CellSystem(sys::Union{AbstractSystem,AtomsVector}) = CellSystem(sys, :)
CellSystem(sys::CellSystem) = sys

Base.getindex(sys::CellSystem, i::Int) = sys.base_system[i]
Base.getindex(sys::CellSystem, c::Colon) = sys.base_system[c]

function Base.getindex(sys::CellSystem, x::Symbol)
    if x == :cell_vectors
        cell_vectors(cell(sys))
    elseif x == :periodicity
        return periodicity(cell(sys))
    else
        return sys.base_system[x]
    end
end

Base.keys(::CellSystem) = (:cell_vectors, :periodicity)
Base.length(sys::CellSystem) = length(sys.base_system)


AtomsBase.cell(sys::CellSystem) = sys.cell


function Base.append!(sys1::T, sys2::T) where{T<:CellSystem}
    @argcheck cell(sys1) == cell(sys2)
    Base.append!(sys1.base_system, sys2.base_system)
    return sys1
end

function Base.:+(sys1::T, sys2::T) where{T<:CellSystem}
    @argcheck cell(sys1) == cell(sys2)
    tmp = deepcopy(sys1)
    append!(tmp, sys2)
    return tmp
end