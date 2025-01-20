

mutable struct SimpleSystem{D, LU, TP} <: AbstractSimpleSystem{D, LU}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    function SimpleSystem(species::AbstractVector{ChemicalSpecies}, r::AbstractVector{<:AbstractVector})
        @argcheck length(species) == length(r)
        D  = (length∘eltype)(r)
        TP = (eltype∘eltype)(r)
        LU = unit(TP)
        new{D, LU, TP}(species, r)
    end
end

function SimpleSystem(species::ChemicalSpecies, pos::AbstractVector{<:Unitful.Length})
    return SimpleSystem([species], [pos])
end

SimpleSystem(sys::SimpleSystem) = sys

function SimpleSystem(ss::Union{AbstractSystem,AtomsVector}, i)
    return SimpleSystem(species(ss, i), position(ss, i))
end

SimpleSystem(sys::Union{AbstractSystem, AtomsVector}) = SimpleSystem(sys, :)

function SimpleSystem(ss::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies)
    i = species(ss, :) .== spc
    return SimpleSystem(ss, i)
end


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

AtomsBase.atomkeys(::SimpleSystem) = (:position, :species)
AtomsBase.hasatomkey(ss::AbstractSimpleSystem, x::Symbol) = x in atomkeys(ss)
AtomsBase.cell(::SimpleSystem{D, LU, TP}) where{D, LU, TP} = IsolatedCell(D, TP)


function Base.append!(sys1::SimpleSystem{D, LU, TP}, sys2::SimpleSystem{D, LU, TP}) where{D, LU, TP}
    Base.append!(sys1.position, sys2.position)
    Base.append!(sys1.species, sys2.species)
    return sys1
end

function Base.deleteat!(sys::SimpleSystem, i)
    Base.deleteat!(sys.species, i)
    Base.deleteat!(sys.position, i)
    return sys
end


##

mutable struct SimpleVelocitySystem{D, LU, UV, TP, TV} <: AbstractSimpleSystem{D, LU}
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
        LU = unit(TP)

        TV = (eltype∘eltype)(v)
        UV = unit(TV)
        new{D, LU, UV, TP, TV}(species, r, v)
    end
end


function SimpleVelocitySystem(sys::Union{AbstractSystem, AtomsVector}, i)
    if hasatomkey(sys, :velocity)
        scp = species(sys, i)
        pos = position(sys, i)
        vel = velocity(sys, i)
        return SimpleVelocitySystem(scp, pos, vel)
    else
        return SimpleSystem(sys, i)
    end
end

SimpleVelocitySystem(sys::Union{AbstractSystem, AtomsVector}) = SimpleVelocitySystem(sys, :)

function SimpleVelocitySystem(
    species::ChemicalSpecies, 
    pos::AbstractVector{<:Unitful.Length}, 
    vel::AbstractVector{<:Unitful.Velocity}
)
    return SimpleVelocitySystem([species], [pos], [vel])
end

function SimpleVelocitySystem(ss::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies)
    i = species(ss, :) .== spc
    return SimpleVelocitySystem(ss, i)
end

Base.getindex(ss::SimpleVelocitySystem, i::Int) = SimpleAtom(ss.species[i], ss.position[i]; velocity=ss.velocity[i])

AtomsBase.atomkeys(::SimpleVelocitySystem) = (:position, :velocity, :species)
AtomsBase.cell(::SimpleVelocitySystem{D, LU, UV, TP, TV}) where{D, LU,UV, TP, TV} = IsolatedCell(D, TP)
AtomsBase.velocity(sys::SimpleVelocitySystem, i) = sys.velocity[i]

function Base.append!(
    sys1::SimpleVelocitySystem{D, LU, UV, TP, TV}, 
    sys2::SimpleVelocitySystem{D, LU, UV, TP, TV}
) where{D, LU, UV, TP, TV}
    Base.append!(sys1.position, sys2.position)
    Base.append!(sys1.species, sys2.species)
    Base.append!(sys1.velocity, sys2.velocity)
    return sys1
end


function Base.deleteat!(sys::SimpleVelocitySystem, i)
    Base.deleteat!(sys.species, i)
    Base.deleteat!(sys.position, i)
    Base.deleteat!(sys.velocity, i)
    return sys
end