


mutable struct AtomicPropertySystem{D, UL, TB, TP} <: AbstractIsolatedSystem{D, UL}
    base_system::TB
    atom_properties::Vector{TP}
    function AtomicPropertySystem(
        sys::AbstractSimpleSystem{D, UL}, 
        properties::AbstractVector{<:NamedTuple}
    ) where {D, UL}
        @argcheck length(sys) == length(properties)
        new{D, UL, typeof(sys), eltype(properties)}(sys, properties)
    end
end

function AtomicPropertySystem(sys::AbstractSystem, properties::NamedTuple)
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

function AtomicPropertySystem(sys::AbstractSystem)
    prop_names = collect( atomkeys(sys) ) # Need Vector for filter to work
    no_special_mass = all( 1:length(sys) ) do i
        mass(sys, i) ≈ mass(species(sys, i))
    end
    if no_special_mass
        filter!(x-> x != :mass, prop_names)
    end
    filter!(x-> ! in(x, (:species, :position, :velocity)), prop_names)
    prop_names = Tuple(prop_names)
    base_sys = hasatomkey(sys, :velocity) ? SimpleVelocitySystem(sys) : SimpleSystem(sys)
    if length(prop_names) == 0
        return base_sys
    end
    a = sys[1]
    el_types = Tuple( typeof(a[key]) for key in prop_names )
    TP = NamedTuple{ prop_names, Tuple{el_types...} }
    prop = Vector{TP}(undef, length(sys))
    for i in 1:length(sys)
        a = sys[i]
        prop[i] = NamedTuple( key =>  a[key] for key in prop_names )
    end
    return AtomicPropertySystem(base_sys, prop)
end


function Base.getindex(sys::AtomicPropertySystem, i::Int)
    Atom(sys.base_system[i]; sys.atom_properties[i]...)
end

Base.getindex(sys::AtomicPropertySystem, x::Symbol) = sys.base_system[x]
Base.keys(::AtomicPropertySystem) = (:cell_vectors, :periodicity)
Base.length(sys::AtomicPropertySystem) = length(sys.base_system)

function AtomsBase.atomkeys(sys::AtomicPropertySystem)
    base_keys = atomkeys(sys.base_system)
    property_keys = keys(sys.atom_properties[1])
    if :mass in property_keys
        base_keys = Tuple( x for x in base_keys if x !=:mass  )
    end
    return (base_keys..., property_keys...)
end

AtomsBase.cell(sys::AtomicPropertySystem) = cell(sys.base_system)
