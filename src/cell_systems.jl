
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

##

function AtomsBase.set_cell_vectors!(sys::CellSystem{D}, cell_matrix::SMatrix{D,D, <:Unitful.Length}) where{D}
    fpos = fractional_coordinates_as_matrix(sys, :)
    new_pos = cell_matrix * fpos
    AtomsBase.set_position!(sys, :, new_pos)
    sys.cell = PeriodicCell( Tuple( x for x in eachcol(cell_matrix) ), periodicity(sys) )
    return sys
end

function AtomsBase.set_cell_vectors!(sys::CellSystem{D}, new_vectors::Vararg{ AbstractVector{<:Unitful.Length} , D}) where{D}
    nc = reduce(hcat, new_vectors)
    return AtomsBase.set_cell_vectors!(sys, nc)
end

function AtomsBase.set_periodicity!(sys::CellSystem{D}, periodicity::NTuple{D, Bool}) where{D}
    sys.cell = PeriodicCell( cell_vectors(sys), periodicity )
    return sys
end

function AtomsBase.set_cell!(sys::CellSystem{D}, ::IsolatedCell{D}) where{D}
    sys = sys.base_system
    return sys
end

function AtomsBase.set_cell!(sys::CellSystem{D}, cell::PeriodicCell{D}) where{D}
    AtomsBase.set_cell_vectors!(sys, cell_matrix(cell))
    return sys
end