
mutable struct CellSystem{D, UL, TB, TC} <: AbstractSystem{D}
    base_system::TB
    cell::TC
    function CellSystem(sys::AbstractIsolatedSystem{D, UL}, cell::PeriodicCell{D,T}) where {D, UL, T}
        new{D, UL, typeof(sys), typeof(cell)}(sys, cell)
    end
end


function CellSystem(sys::AbstractSystem)
    base_sys = AtomicPropertySystem(sys)
    if isa(cell(sys), IsolatedCell)
        return base_sys
    end
    return CellSystem(base_sys, cell(sys))
end

Base.getindex(cs::CellSystem, i::Int) = cs.base_system[i]
Base.getindex(cs::CellSystem, c::Colon) = cs.base_system[c]
Base.getindex(cs::CellSystem, x::Symbol) = cs.base_system[x]
Base.keys(::CellSystem) = (:cell_vectors, :periodicity)
Base.length(cs::CellSystem) = length(cs.base_system)

AtomsBase.cell(cs::CellSystem) = cs.cell