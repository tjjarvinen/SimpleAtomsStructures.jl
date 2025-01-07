


mutable struct AtomicPropertySystem{D, LU, TB, TP, SM} <: AbstractIsolatedSystem{D, LU}
    base_system::TB
    atom_properties::Vector{TP}
    function AtomicPropertySystem(
        sys::AbstractSimpleSystem{D, LU}, 
        properties::AbstractVector{<:NamedTuple}
    ) where {D, LU}
        @argcheck length(sys) == length(properties)
        new{D, LU, typeof(sys), eltype(properties), haskey(properties[1], :mass)}(sys, properties)
    end
end

function AtomicPropertySystem(sys::AbstractSystem, properties::NamedTuple)
    # ignore atom_properties on sys and use properties as atom_properties
    # also ignore global system properties
    @argcheck all(x->length(x)==length(sys), properties)
    prop_names = Tuple( x  for x in keys(properties) if !(x in (:species, :position, :velocity)) )
    el_types = Tuple( eltype(properties[key]) for key in prop_names )
    TP = NamedTuple{ prop_names, Tuple{el_types...} }
    prop = Vector{TP}(undef, length(sys))
    for i in 1:length(sys)
        prop[i] = NamedTuple( key =>  properties[key][i] for key in prop_names )
    end
    return AtomicPropertySystem(sys, prop)
end


AtomicPropertySystem(sys::AtomicPropertySystem) = sys
function AtomicPropertySystem(sys::AtomicPropertySystem, i)
    tmp = SimpleVelocitySystem(sys.base_system, i)
    return AtomicPropertySystem( tmp, sys.atom_properties[i] )
end

function AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, i)
    prop_names = [ k for k in atomkeys(sys) if ! in(k, (:species, :position, :velocity)) ]
    if length(prop_names) == 0
        return SimpleVelocitySystem(sys, i)
    end
    tmp = map( sys[i] ) do at
        NamedTuple( k=>at[k] for k in prop_names )
    end
    return AtomicPropertySystem( SimpleVelocitySystem(sys, i), tmp )
end

AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, ::Colon) = AtomicPropertySystem(sys, 1:length(sys))
AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector})          = AtomicPropertySystem(sys, :)
AtomicPropertySystem(sys::AbstractSimpleSystem, i) = SimpleVelocitySystem(sys, i)
AtomicPropertySystem(sys::AbstractSimpleSystem)    = SimpleVelocitySystem(sys)

function AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, i::BitVector)
    @argcheck length(sys) == length(i)
    j = (1:length(sys))[i]
    return AtomicPropertySystem(sys, j)
end

function AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, spc::ChemicalSpecies)
    i = species(sys, :) .== spc
    return AtomicPropertySystem(sys, i)
end


function Base.getindex(sys::AtomicPropertySystem, i::Int)
    tmp = sys.base_system[i]
    SimpleAtom( (; tmp.data..., sys.atom_properties[i]...) )
end

Base.getindex(sys::AtomicPropertySystem, x::Symbol) = sys.base_system[x]
Base.keys(::AtomicPropertySystem) = (:cell_vectors, :periodicity)
Base.length(sys::AtomicPropertySystem) = length(sys.base_system)

function AtomsBase.atomkeys(sys::AtomicPropertySystem)
    base_keys = AtomsBase.atomkeys(sys.base_system)
    property_keys = _property_keys(sys)
    # remove double mass
    if :mass in property_keys
        base_keys = Tuple( x for x in base_keys if x !=:mass  )
    end
    return (base_keys..., property_keys...)
end

_property_keys(sys::AtomicPropertySystem) = keys(sys.atom_properties[1])
_has_property_key(sys::AtomicPropertySystem, x::Symbol) = in(x, _property_keys(sys)) 


function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, i::Int) where {D, LU, TB, TP}
    return getindex(sys.atom_properties[i], :mass)
end
function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, i) where {D, LU, TB, TP}
    return [ mass(sys, j)  for j in i ]
end
function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, ::Colon) where {D, LU, TB, TP}
    return mass(sys, 1:length(sys))
end

function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, false}, i) where {D, LU, TB, TP}
    return mass(sys.base_system, i)
end

function Base.append!(sys1::T, sys2::T) where{T<:AtomicPropertySystem}
    Base.append!(sys1.base_system, sys2.base_system)
    Base.append!(sys1.atom_properties, sys2.atom_properties)
    return sys1
end