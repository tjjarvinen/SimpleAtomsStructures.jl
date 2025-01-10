
mutable struct GeneralSystem{D, LU, TB} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    properties::Dict{Symbol, Any}
    function GeneralSystem(
        sys::Union{CellSystem{D,LU,TT,TC}, AbstractIsolatedSystem{D,LU}};
        kwargs...
    ) where {D,LU,TT,TC}
        props = Dict{Symbol, Any}( k=>v for (k,v) in pairs(kwargs) if ! in(k, (:cell, :cell_vectors, :periodicity)) )
        new{D, LU, typeof(sys)}(sys, props)
    end
end

GeneralSystem(sys::GeneralSystem) = sys

Base.getindex(sys::GeneralSystem, i::Int) = sys.base_system[i]
Base.getindex(sys::GeneralSystem, c::Colon) = sys.base_system[c]

function Base.getindex(sys::GeneralSystem, x::Symbol)
    if x in (:cell_vectors, :periodicity)
        return sys.base_system[x]
    else
        return sys.properties[x]
    end
end

function Base.keys(sys::GeneralSystem)
    tmp = keys(sys.properties)
    return (:cell_vectors, :periodicity, tmp...)
end

Base.haskey(sys::GeneralSystem, x::Symbol) = in(x, keys(sys))


## Generic builders

function GenericSystem(sys::AbstractSystem; kwargs...)
    # Figure out do we need to update the cell
    new_cell = nothing
    if haskey(kwargs, :cell_vectors) && haskey(kwargs, :periodicity)
        new_cell = PeriodicCell( kwargs[:cell_vectors] )
    elseif haskey(kwargs, :cell_vectors) && isa(cell(sys), PeriodicCell)
        new_cell = PeriodicCell( kwargs[:cell_vectors], periodicity(sys) )
    elseif haskey(kwargs, :periodicity) && isa(cell(sys), PeriodicCell)
        new_cell = PeriodicCell( cell_vectors(sys), kwargs[:periodicity] )
    elseif isa(cell(sys), IsolatedCell) &&
            ( haskey(kwargs, :cell_vectors) || haskey(kwargs, :periodicity))
        throw( ArgumentError("Not enough information to update cell") )
    end

    if isnothing(new_cell)
        tsys = CellSystem(sys)
        if length(kwargs) == 0 && length(keys(sys)) == 2 # only :cell_vectors and :periodicity
            return tsys
        elseif length(kwargs) == 0 && length(keys(sys)) > 2
            tmp = NamedTuple( k=>sys[k] for k in keys(sys) )
            return GeneralSystem(tsys; tmp...)
        elseif length(kwargs) > 0 && length(keys(sys)) == 2
            return GeneralSystem(tsys; kwargs...)
        else
            tmp = Dict{Symbol, Any}( k=>sys[k] for k in keys(sys) )
            foreach( pairs(kwargs) ) do (k,v)
                tmp[k] = v
            end
            return GeneralSystem(tsys; tmp...)
        end
    end

    tsys = CellSystem( AtomicPropertySystem(sys), new_cell )
    if length(kwargs) > 0 && length(keys(sys)) == 2
        return GeneralSystem(tsys; kwargs...)
    else
        tmp = Dict{Symbol, Any}( k=>sys[:k] for k in keys(sys) )
        foreach( pairs(kwargs) ) do (k,v)
            tmp[k] = v
        end
        return GeneralSystem(tsys; tmp...)
    end
end

function GenericSystem(sys::AtomsVector; kwargs...)
    tmp = AtomicPropertySystem(sys)
    if length(kwargs) > 0
        return GenericSystem(tmp; kwargs...)
    end
    return tmp
end

function GenericSystem(sys::AbstractVector{<:Union{AtomsBase.Atom, AtomsBase.AtomView}}; kwargs...)
    tmp = SimpleAtom.(sys)
    return GenericSystem(tmp; kwargs...)
end

function GenericSystem(
    species::AbstractVector{<:AtomsBase.AtomId},
    r::AbstractVector{<:AbstractVector};
    kwargs...
)   
    spc = ChemicalSpecies.(species)
    tmp = SimpleSystem(spc, r)
    if length(kwargs) > 0
        return GenericSystem(tmp; kwargs...)
    end
    return tmp
end

function GenericSystem(
    species::AbstractVector{<:AtomsBase.AtomId}, 
    r::AbstractVector{<:AbstractVector},
    v::AbstractVector{<:AbstractVector}; 
    kwargs...
)   
    spc = ChemicalSpecies.(species)
    tmp = SimpleVelocitySystem(spc, r, v)
    if length(kwargs) > 0
        return GenericSystem(tmp; kwargs...)
    end
    return tmp
end


GenericSystem(sys::AbstractSystem, i) = CellSystem(sys, i)
