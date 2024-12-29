
struct SimpleAtom{D, TD}
    data::TD
    function SimpleAtom(
        spc::ChemicalSpecies, 
        r::AbstractVector{<:Unitful.Length}; 
        kwords...
    )   
        if haskey(kwords, :velocity) && length(kwords[:velocity]) != length(r)
            throw( ArgumentError("Position and velocity have different dimensions") )
        end
        tmp = ( species=spc, position=SVector(r...), kwords... )
        new{length(r), typeof(tmp)}(tmp)
    end
    function SimpleAtom(data::NamedTuple)
        @argcheck haskey(data, :species)
        @argcheck haskey(data, :position)
        @argcheck ! ( haskey(data, :velocity) && length(data.velocity) != length(data.position) )
        new{length(data.position), typeof(data)}(data)
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


AtomsBase.n_dimensions(::SimpleAtom{D, TD}) where {D, TD} = D

AtomsBase.velocity(sa::SimpleAtom) = haskey(sa, :mass) ? sa.data.velocity : missing
AtomsBase.position(sa::SimpleAtom) = sa.data.position
AtomsBase.mass(sa::SimpleAtom)     = haskey(sa, :mass) ? sa.data.mass : AtomsBase.mass(AtomsBase.species(sa))
AtomsBase.species(sa::SimpleAtom)  = sa.data.species

AtomsBase.atom_name(sa::SimpleAtom)     = atom_name(species(sa))
AtomsBase.atomic_symbol(sa::SimpleAtom) = atomic_symbol(species(sa))
AtomsBase.atomic_number(sa::SimpleAtom) = atomic_number(species(sa))
AtomsBase.element(sa::SimpleAtom)       = element(species(sa))

Base.getindex(sa::SimpleAtom, x::Symbol) = sa.data[x]
Base.haskey(sa::SimpleAtom, x::Symbol)   = Base.haskey(sa.data, x)
Base.pairs(sa::SimpleAtom)               = pairs(sa.data)

Base.show(io::IO, sa::SimpleAtom) = AtomsBase.show_atom(io, sa)
Base.show(io::IO, mime::MIME"text/plain", sa::SimpleAtom) = AtomsBase.show_atom(io, mime, sa)