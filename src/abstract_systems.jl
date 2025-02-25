
abstract type AbstractCompositeSystem{D, LU} <: AtomsBase.AbstractSystem{D}    end
abstract type AbstractIsolatedSystem{D, LU}  <: AbstractCompositeSystem{D, LU} end
abstract type AbstractSimpleSystem{D, LU}    <: AbstractIsolatedSystem{D, LU}  end

abstract type AbstractTrajectory{D, LU} <: AbstractVector{AbstractSystem{D}} end

## Defive some properties

AtomsBase.atomkeys(sys::AbstractCompositeSystem) = AtomsBase.atomkeys(sys.base_system)
AtomsBase.cell(sys::AbstractCompositeSystem) = AtomsBase.cell(sys.base_system)
AtomsBase.mass(sys::AbstractCompositeSystem, i) = AtomsBase.mass(sys.base_system, i)
AtomsBase.position(sys::AbstractCompositeSystem, i) =  position(sys.base_system, i)
AtomsBase.species(sys::AbstractCompositeSystem, i) = AtomsBase.species(sys.base_system, i)
AtomsBase.velocity(sys::AbstractCompositeSystem, i) =  AtomsBase.velocity(sys.base_system, i)

AtomsBase.set_position!(sys::AbstractCompositeSystem, i, x) = AtomsBase.set_position!(sys.base_system, i, x)

Base.length(sys::AbstractCompositeSystem) = length(sys.base_system)

function Base.:+(sys1::T, sys2::T) where{T<:AbstractIsolatedSystem}
    tmp = deepcopy(sys1)
    append!(tmp, sys2)
    return tmp
end


AtomsBase.mass(ss::AbstractSimpleSystem, i) = mass.(ss.species[i])
AtomsBase.position(ss::AbstractSimpleSystem, i::Int) = ss.position[i]
AtomsBase.position(ss::AbstractSimpleSystem, i) = view(ss.position, i)
AtomsBase.position(ss::AbstractSimpleSystem, ::Colon) = ss.position
AtomsBase.species(ss::AbstractSimpleSystem, i::Int) = ss.species[i]
AtomsBase.species(ss::AbstractSimpleSystem, i) = view(ss.species, i)
AtomsBase.species(ss::AbstractSimpleSystem, ::Colon) = ss.species
AtomsBase.element_symbol(ss::AbstractSimpleSystem, i) = element_symbol.(species(ss, i))

function AtomsBase.set_species!(sys::AbstractSimpleSystem, i, x)
    setindex!(sys.species, x, i)
    return sys
end

function AtomsBase.set_position!(sys::AbstractSimpleSystem, i, x)
    setindex!(sys.position, x, i)
    return sys
end

function AtomsBase.set_position!(sys::AbstractSimpleSystem{D}, ::Colon, x::AbstractMatrix{<:Unitful.Length}) where{D}
    @argcheck size(x, 2) == length(sys)
    tmp = reinterpret(reshape, SVector{D, eltype(x)}, x)
    AtomsBase.set_position!(sys, :, tmp)
    return sys
end