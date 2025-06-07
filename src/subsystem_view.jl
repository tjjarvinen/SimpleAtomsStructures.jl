
mutable struct SimpleSystemView{D, LU, TP, TI} <: AbstractSimpleSystem{D, LU}
    species::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, true}
    position::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, true}
    function SimpleSystemView(sys::AbstractCompositeSystem{D, LU}, i) where {D, LU}
        vspc = view( species(sys, :), i)
        vpos = view( position(sys, :), i)
        T = typeof(vspc)
        TI = T.parameters[4]
        TP = (eltype ∘ eltype)(vpos) 
        new{D, LU, TP, TI}(vspc, vpos)
    end
end


mutable struct SimpleVelocitySystemView{D, LU, TP, TV, TI} <: AbstractSimpleSystem{D, LU}
    species::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, true}
    position::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, true}
    velocity::SubArray{SVector{D,TV}, 1, Vector{SVector{D,TV}}, TI, true}
    function SimpleVelocitySystemView(sys::AbstractCompositeSystem{D, LU}, i) where {D, LU}
        vspc = view( species(sys, :), i)
        vpos = view( position(sys, :), i)
        vvel = view( velocity(sys, :), i)
        T = typeof(vspc)
        TI = T.parameters[4]
        TP = (eltype ∘ eltype)(vpos)
        TV = (eltype ∘ eltype)(vvel)    
        new{D, LU, TP, TV, TI}(vspc, vpos, vvel)
    end
end

mutable struct AtomicPropertySystemView{D, LU, TB, TP, TI, SM} <: AbstractIsolatedSystem{D, LU}
    base_system::TB
    atom_properties::SubArray{TP, 1, Vector{TP}, TI, true}
    function AtomicPropertySystemView(sys::AtomicPropertySystem{D, LU, TT, TP, SM}, i) where {D, LU, TT, TP, SM}
        base = system_view(sys.base_system, i)
        atom_properties = view(sys.atom_properties, i)
        T = typeof(atom_properties)
        TI = T.parameters[4]
        new{D, LU, typeof(base), TP, TI, SM}(base, atom_properties)
    end
end

mutable struct CellSystemView{D, LU, TB, TC} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    cell::TC
    function CellSystemView(sys::CellSystem{D, LU, TT, TC}, i) where {D, LU, TT, TC}
        base = system_view(sys.base_system, i)
        new{D, LU, typeof(base), TC}(base, sys.cell)
    end
end

##

Base.getindex(sys::SimpleSystemView, i::Int) = SimpleAtom( sys.species[i],  sys.position[i] )
Base.getindex(sys::SimpleVelocitySystemView, i::Int) = SimpleAtom( sys.species[i],  sys.position[i], sys.velocity[i] )
function Base.getindex(sys::AtomicPropertySystemView, i::Int)
    tmp = sys.base_system[i]
    SimpleAtom( (; tmp.data..., sys.atom_properties[i]...) )
end

function Base.getindex(sys::CellSystemView, x::Symbol)
    if x == :cell_vectors
        cell_vectors(cell(sys))
    elseif x == :periodicity
        return periodicity(cell(sys))
    else
        return sys.base_system[x]
    end
end

AtomsBase.atomkeys(::SimpleSystemView) = (:position, :species)
AtomsBase.atomkeys(::SimpleVelocitySystemView) = (:position, :species, :velocity)

AtomsBase.cell(::SimpleSystemView{D, LU, TP}) where {D, LU, TP} = IsolatedCell(D, TP)
AtomsBase.cell(::SimpleVelocitySystemView{D, LU, TP}) where {D, LU, TP} = IsolatedCell(D, TP)
AtomsBase.cell(sys::CellSystemView) = sys.cell


system_view(sys::Union{SimpleSystem, SimpleSystemView}, i) = SimpleSystemView(sys, i)
system_view(sys::Union{SimpleVelocitySystem, SimpleVelocitySystemView}, i) = SimpleVelocitySystemView(sys, i)
system_view(sys::Union{AtomicPropertySystem, AtomicPropertySystemView}, i) = AtomicPropertySystemView(sys, i)
system_view(sys::Union{CellSystem, CellSystemView}, i) = CellSystemView(sys, i)

system_view(sys::GeneralSystem, i) = CellSystemView(sys.base_system, i)  
    
