
# This might not be the most efficient implementation due to
# check that position and velocity have same dimensions
# so it might need an update
struct SimpleAtom{D, TD, TP}
    data::TD
    function SimpleAtom(
        spc::ChemicalSpecies, 
        r::AbstractVector{<:Unitful.Length}; 
        kwargs...
    )   
        if haskey(kwargs, :velocity) && length(kwargs[:velocity]) != length(r)
            throw( ArgumentError("Position and velocity have different dimensions") )
        end
        tmp = ( species=spc, position=SVector(r...), kwargs... )
        new{length(r), typeof(tmp), eltype(r)}(tmp)
    end
    function SimpleAtom(data::NamedTuple)
        @argcheck haskey(data, :species)
        @argcheck haskey(data, :position)
        @argcheck ! ( haskey(data, :velocity) && length(data.velocity) != length(data.position) )
        new{length(data.position), typeof(data), eltype(data.position)}(data)
    end
end

function SimpleAtom(at::Union{AtomsBase.Atom, AtomsBase.AtomView})
    spc = species(at)
    r = position(at)
    properties = NamedTuple( k=>at[k] for k in keys(at) if ! in(k, (:species, :position)))
    return SimpleAtom(spc, r; properties...)
end

function SimpleAtom(sa::SimpleAtom; kwargs...)
    length(kwargs) == 0 && return sa
    tmp = Dict{Symbol, Any}( pairs(sa) )
    foreach( pairs(kwargs) ) do (k,v)
        tmp[k] = v
    end
    return SimpleAtom( NamedTuple( tmp ) )
end

function SimpleAtom(id::AtomsBase.AtomId, r::AbstractVector{<:Unitful.Length}; kwargs...)
    return SimpleAtom(ChemicalSpecies(id), r; kwargs...)
end

function SimpleAtom(pr::Pair; kwargs...)
    return SimpleAtom(pr[1], pr[2]; kwargs...)
end

AtomsBase.n_dimensions(::SimpleAtom{D, TD}) where {D, TD} = D

AtomsBase.velocity(sa::SimpleAtom) = haskey(sa, :velocity) ? sa.data.velocity : missing
AtomsBase.position(sa::SimpleAtom) = sa.data.position
AtomsBase.mass(sa::SimpleAtom)     = haskey(sa, :mass) ? sa.data.mass : AtomsBase.mass(AtomsBase.species(sa))
AtomsBase.species(sa::SimpleAtom)  = sa.data.species

AtomsBase.atom_name(sa::SimpleAtom)     = atom_name(species(sa))
AtomsBase.atomic_symbol(sa::SimpleAtom) = atomic_symbol(species(sa))
AtomsBase.atomic_number(sa::SimpleAtom) = atomic_number(species(sa))
AtomsBase.element(sa::SimpleAtom)       = element(species(sa))

Base.getindex(sa::SimpleAtom, x::Symbol) = x == :mass ? mass(sa) : sa.data[x]
Base.haskey(sa::SimpleAtom, x::Symbol)   = Base.haskey(sa.data, x)
Base.keys(sa::SimpleAtom)                = Base.keys(sa.data)
Base.pairs(sa::SimpleAtom)               = pairs(sa.data)

Base.show(io::IO, sa::SimpleAtom) = AtomsBase.show_atom(io, sa)
Base.show(io::IO, mime::MIME"text/plain", sa::SimpleAtom) = AtomsBase.show_atom(io, mime, sa)

Base.convert(SimpleAtom, at::Union{AtomsBase.Atom, AtomsBase.AtomView}) = SimpleAtom(at)

## Vector properties

const AtomsVector{D,TD, TP} = AbstractVector{SimpleAtom{D, TD, TP}}

AtomsBase.atomkeys(sys::AtomsVector) = keys(sys[1])
AtomsBase.hasatomkey(sys::AtomsVector, x::Symbol) = haskey(sys[1], x)

AtomsBase.cell(::AtomsVector{D, TD, TP}) where{D, TD, TP} = IsolatedCell(D, TP)
AtomsBase.cell_vectors(sys::AtomsVector) = cell_vectors(cell(sys))
AtomsBase.periodicity(sys::AtomsVector) = periodicity(cell(sys))

AtomsBase.mass(sys::AtomsVector, i) = map(j->sys[j][:mass], i)
AtomsBase.mass(sys::AtomsVector, ::Colon) = map(j->sys[j][:mass], eachindex(sys))

AtomsBase.position(sys::AtomsVector, i) = map(j->sys[j][:position], i)
AtomsBase.position(sys::AtomsVector, ::Colon) = map(j->sys[j][:position], eachindex(sys))

AtomsBase.species(sys::AtomsVector, i) = map(j->sys[j][:species], i)
AtomsBase.species(sys::AtomsVector, ::Colon) = map(j->sys[j][:species], eachindex(sys))

AtomsBase.velocity(sys::AtomsVector, i) = map(j->sys[j][:velocity], i)
AtomsBase.velocity(sys::AtomsVector, ::Colon) = map(j->sys[j][:velocity], eachindex(sys))
