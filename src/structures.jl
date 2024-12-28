
abstract type AbstractSimpleSystem{D, UL, TF}  <: AtomsBase.AbstractSystem{D} end

mutable struct SimpleSystem{D, UL, TF, TP} <: AbstractSimpleSystem{D, UL, TF}
    species::Vector{ChemicalSpecies}
    positions::Vector{SVector{D, TP}}
    function SimpleSystem(species::AbstractVector{ChemicalSpecies}, r::AbstractVector{<:AbstractVector})
        @argcheck length(species) == length(r)
        D = (length∘eltype)(r)
        TP = (eltype∘eltype)(r)
        TF = (typeof∘ustrip)(zero(TP))
        UL = unit(TP)
        new{D, UL, TF, TP}(species, r)
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


Base.getindex(ss::SimpleSystem, i::Int) = Atom(ss.species[i], ss.positions[i])
function Base.getindex(ss::SimpleSystem, x::Symbol)
    if x === :cell_vectors
        return cell_vectors(ss)
    elseif x === :periodicity
        return periodicity(ss)
    else
        throw( KeyError(x) )
    end
end


Base.keys(::SimpleSystem) = (:cell_vectors, :periodicity)
Base.length(ss::SimpleSystem) = length(ss.species)

AtomsBase.atomkeys(::SimpleSystem) = (:position, :species, :mass)
AtomsBase.cell(::SimpleSystem{D, UL, TF, TP}) where{D, UL, TF, TP} = IsolatedCell(D, TP)
AtomsBase.mass(ss::SimpleSystem, i) = mass.(ss.species[i])
AtomsBase.position(ss::SimpleSystem, i) = ss.positions[i]
AtomsBase.species(ss::SimpleSystem, i) = ss.species[i]
AtomsBase.element_symbol(ss::SimpleSystem, i) = element_symbol.(species(ss, i))


function get_subsystem(ss::SimpleSystem, i)
    return SimpleSystem(ss.species[i], ss.positions[i])
end

function get_subsystem(ss::SimpleSystem, spc::ChemicalSpecies)
    i = ss.species .== spc
    return SimpleSystem(ss.species[i], ss.positions[i])
end

##

# mutable struct SimpleVelocitySystem{D, UL, UT, TF, TP, TV} <: AbstractSimpleSystem{D}
#     base_system::SimpleSystem{D, UL, TF, TP}
#     velocities::Vector{SVector{D, TV}}
#     function SimpleVelocitySystem(
#         sys::SimpleSystem{D, UL, TF, TP},
#         v::AbstractVector{<:AbstractVector}
#     ) where {D, UL, TF, TP}
#         @argcheck length(sys) == length(v)
#         @argcheck D == (length∘eltype)(v)
#         _TV = (eltype∘eltype)(v)
#         UV = unit(_TV)
#         UT = UL / UV |> upreferred
#         TV = typeof( zero(TF) * UL/UT )
#         new{D, UL, UT, TF, TP, TV}()(sys, v)
#     end
# end