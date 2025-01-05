
mutable struct CellSystem{D, LU, TB, TC} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    cell::TC
    function CellSystem(sys::AbstractIsolatedSystem{D, LU}, cell::PeriodicCell{D,T}) where {D, LU, T}
        new{D, LU, typeof(sys), typeof(cell)}(sys, cell)
    end
end


function CellSystem(sys::AbstractSystem)
    base_sys = AtomicPropertySystem(sys)
    if isa(cell(sys), IsolatedCell)
        return base_sys
    end
    return CellSystem(base_sys, cell(sys))
end

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