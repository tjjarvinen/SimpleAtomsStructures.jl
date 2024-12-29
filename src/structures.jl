
abstract type AbstractIsolatedSystem{D, UL}  <: AtomsBase.AbstractSystem{D} end
abstract type AbstractSimpleSystem{D, UL}  <: AbstractIsolatedSystem{D, UL} end

mutable struct SimpleSystem{D, UL, TP} <: AbstractSimpleSystem{D, UL}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    function SimpleSystem(species::AbstractVector{ChemicalSpecies}, r::AbstractVector{<:AbstractVector})
        @argcheck length(species) == length(r)
        D  = (length∘eltype)(r)
        TP = (eltype∘eltype)(r)
        UL = unit(TP)
        new{D, UL, TP}(species, r)
    end
end

function SimpleSystem(sys::AbstractSystem)
    scp = species(sys, :)
    pos = position(sys, :)
    return SimpleSystem(scp, pos)
end

function SimpleSystem(species::ChemicalSpecies, pos::AbstractVector{<:Unitful.Length})
    return SimpleSystem([species], [pos])
end

SimpleSystem(sys::SimpleSystem) = sys


Base.getindex(ss::SimpleSystem, i::Int) = SimpleAtom(ss.species[i], ss.position[i])
function Base.getindex(ss::AbstractSimpleSystem, x::Symbol)
    if x === :cell_vectors
        return cell_vectors(ss)
    elseif x === :periodicity
        return periodicity(ss)
    else
        throw( KeyError(x) )
    end
end


Base.haskey(ss::AbstractSimpleSystem, x::Symbol) = in(x, keys(ss) )
Base.keys(::AbstractSimpleSystem) = (:cell_vectors, :periodicity)
Base.length(ss::AbstractSimpleSystem) = length(ss.species)

AtomsBase.atomkeys(::SimpleSystem) = (:position, :species, :mass)
AtomsBase.hasatomkey(ss::AbstractSimpleSystem, x::Symbol) = x in atomkeys(ss)
AtomsBase.cell(::SimpleSystem{D, UL, TP}) where{D, UL, TP} = IsolatedCell(D, TP)
AtomsBase.mass(ss::AbstractSimpleSystem, i) = mass.(ss.species[i])
AtomsBase.position(ss::AbstractSimpleSystem, i) = ss.position[i]
AtomsBase.species(ss::AbstractSimpleSystem, i) = ss.species[i]
AtomsBase.element_symbol(ss::AbstractSimpleSystem, i) = element_symbol.(species(ss, i))


function get_subsystem(ss::SimpleSystem, i)
    return SimpleSystem(ss.species[i], ss.position[i])
end

function get_subsystem(ss::AbstractSimpleSystem, spc::ChemicalSpecies)
    i = ss.species .== spc
    return get_subsystem(ss, i)
end

##

mutable struct SimpleVelocitySystem{D, UL, UV, TP, TV} <: AbstractSimpleSystem{D, UL}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    velocity::Vector{SVector{D, TV}}
    function SimpleVelocitySystem(
        species::AbstractVector{ChemicalSpecies}, 
        r::AbstractVector{<:AbstractVector}, 
        v::AbstractVector{<:AbstractVector}
    )
        @argcheck length(species) == length(r) == length(v)
        @argcheck length( eltype(r) ) == length( eltype(v) )
        D  = (length∘eltype)(r)
        TP = (eltype∘eltype)(r)
        UL = unit(TP)

        TV = (eltype∘eltype)(v)
        UV = unit(TV)
        new{D, UL, UV, TP, TV}(species, r, v)
    end
end


function SimpleVelocitySystem(sys::AbstractSystem)
    @argcheck hasatomkey(sys, :velocity)
    scp = species(sys, :)
    pos = position(sys, :)
    vel = velocity(sys, :)
    return SimpleVelocitySystem(scp, pos, vel)
end

function SimpleVelocitySystem(
    species::ChemicalSpecies, 
    pos::AbstractVector{<:Unitful.Length}, 
    vel::AbstractVector{<:Unitful.Velocity}
)
    return SimpleVelocitySystem([species], [pos], [vel])
end

Base.getindex(ss::SimpleVelocitySystem, i::Int) = SimpleAtom(ss.species[i], ss.position[i]; velocity=ss.velocity[i])

AtomsBase.atomkeys(::SimpleVelocitySystem) = (:position, :velocity, :species, :mass)
AtomsBase.cell(::SimpleVelocitySystem{D, UL, UV, TP, TV}) where{D, UL,UV, TP, TV} = IsolatedCell(D, TP)

function get_subsystem(ss::SimpleVelocitySystem, i)
    return SimpleVelocitySystem(ss.species[i], ss.position[i], ss.velocity[i])
end