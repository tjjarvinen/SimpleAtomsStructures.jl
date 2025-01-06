
abstract type AbstractCompositeSystem{D, LU} <: AtomsBase.AbstractSystem{D} end
abstract type AbstractIsolatedSystem{D, LU}  <: AbstractCompositeSystem{D, LU} end
abstract type AbstractSimpleSystem{D, LU}    <: AbstractIsolatedSystem{D, LU} end


## Defive some properties

AtomsBase.atomkeys(sys::AbstractCompositeSystem) = AtomsBase.atomkeys(sys.base_system)
AtomsBase.cell(sys::AbstractCompositeSystem) = AtomsBase.cell(sys.base_system)
AtomsBase.mass(sys::AbstractCompositeSystem, i) = AtomsBase.mass(sys.base_system, i)
AtomsBase.position(sys::AbstractCompositeSystem, i) =  position(sys.base_system, i)
AtomsBase.species(sys::AbstractCompositeSystem, i) = AtomsBase.species(sys.base_system, i)
AtomsBase.velocity(sys::AbstractCompositeSystem, i) =  AtomsBase.velocity(sys.base_system, i)

Base.length(sys::AbstractCompositeSystem) = length(sys.base_system)

function Base.:+(sys1::T, sys2::T) where{T<:AbstractIsolatedSystem}
    tmp = deepcopy(sys1)
    append!(tmp, sys2)
    return tmp
end