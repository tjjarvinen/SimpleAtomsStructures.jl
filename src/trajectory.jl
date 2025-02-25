

struct SystemView{D, LU, TT} <: AtomsBase.AbstractSystem{D}
    parent_trajectory::TT
    frame::Int
    natoms::Int
    function SystemView(trj::AbstractTrajectory{D, LU}, frame::Int) where{D, LU}
        @argcheck 1 <= frame <= length(trj)
        new{D, LU, typeof(trj)}(trj, frame, length(trj.species))
    end
end

function Base.getindex(sv::SystemView, i::Int)
    @argcheck 1 <= i <= sv.natoms
    j = (i-1)*sv.natoms + i
    return SimpleAtom(sv.parent_trajectory.species[j], sv.parent_trajectory.position[j])
end

function Base.getindex(sv::SystemView, ::Colon)
    i = (i-1)*sv.natoms
    j = i + sv.natoms
    return [SimpleAtom(sv.parent_trajectory.species[k], sv.parent_trajectory.position[k]) for k in i:j]
end

function Base.getindex(sv::SystemView, i::UnitRange)
    @argcheck 1 <= first(i) <= last(i) <= sv.natoms
    return [SimpleAtom(sv.parent_trajectory.species[k], sv.parent_trajectory.position[k]) for k in i]
end

function Base.length(sv::SystemView)
    return sv.natoms
end

function Base.getindex(sv::SystemView, x::Symbol)
    if x === :cell_vectors
        return cell_vectors(sv)
    elseif x === :periodicity
        return periodicity(sv)
    else
        throw( KeyError(x) )
    end
    
end


Base.haskey(sv::SystemView, x::Symbol) = in(x, keys(sv) )
Base.keys(::SystemView) = (:cell_vectors, :periodicity)

AtomsBase.atomkeys(::SystemView) = (:species, :position)
AtomsBase.hasatomkey(::SystemView, key) = key in atomkeys(SystemView)
AtomsBase.cell(sv::SystemView) = cell(sv.parent_trajectory, sv.frame)

AtomsBase.species(sv::SystemView, i) = species(sv.parent_trajectory, i)
AtomsBase.mass(sv::SystemView, i) = mass(sv.parent_trajectory, i)

AtomsBase.position(sv::SystemView, i) = position(sv.parent_trajectory, i, sv.frame)

##

mutable struct Trajectory{D, LU, TP} <: AbstractTrajectory{D, LU}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    natoms::Int
    function Trajectory(
        spc::AbstractVector{ChemicalSpecies},
        pos::AbstractVector{SVector{D, TP}}
    ) where {D, TP<:Unitful.Length}
        @argcheck length(spc) == length(pos)
        LU = unit(TP)
        # need copies to not break things - trajectory is mutable
        new{D, LU, TP}(deepcopy(spc), deepcopy(pos), length(spc))
    end    
end

function Trajectory(sys::AbstractSystem)
    spc = species(sys, :) 
    pos = position(sys, :) 
    return Trajectory(spc, pos)
end

function Trajectory(vsys::AbstractVector{<:AbstractSystem})
    tmp = Trajectory(vsys[1])
    append!(tmp, vsys[2:end])
    return tmp
end

function Base.length(trj::Trajectory)
    return length(trj.position)//length(trj.species) |> Int
end

Base.eltype(::Type{Trajectory{D, LU, TP}}) where{D, LU, TP} = SystemView{D, LU, Trajectory{D, LU, TP}}
Base.size(trj::Trajectory) = (length(trj), )

Base.show(io::IO, trj::Trajectory) = print(io, "Trajectory{", eltype(trj), "} with ", length(trj), " frames")

function Base.push!(trj::Trajectory{D}, sys::AbstractSystem{D}; check_species=false) where{D}
    @argcheck trj.natoms == length(sys)
    if check_species
        @argcheck all( i -> species(trj, i) === species(sys, i) for i = 1:trj.natoms )
    end
    append!(trj.position, position(sys, :))
    return trj
end

function Base.append!(trj::Trajectory{D}, addtraj::AbstractVector{<:AbstractSystem{D}}) where{D}
    # NOTE species are not tested
    @argcheck all( trj.natoms == length(sys) for sys in addtraj )
    pos = mapreduce(vcat, addtraj) do frame
        position(frame, :)
    end
    append!(trj.position, pos)
    return trj
end

function Base.getindex(trj::Trajectory, i::Int)
    @argcheck 1 <= i <= length(trj)
    return SystemView(trj, i)
end

Base.haskey(trj::Trajectory, x::Symbol) = in(x, keys(trj) )
Base.keys(::Trajectory) = (:cell_vectors, :periodicity)

AtomsBase.species(trj::Trajectory, i, frame) = species(trj, i)
AtomsBase.species(trj::Trajectory, i) = view(trj.species, i)
AtomsBase.species(trj::Trajectory, ::Colon) = trj.species
AtomsBase.species(tra::Trajectory, i::Int) = tra.species[i]

AtomsBase.mass(trj::Trajectory, i, frame) = mass(trj, i)
AtomsBase.mass(trj::Trajectory, i) = mass.( species(trj, i) )
AtomsBase.cell(::Trajectory{D, LU, TP}, ::Any) where{D,LU,TP} = IsolatedCell(D, TP)
AtomsBase.periodicity(trj::Trajectory, frame) = periodicity(cell(trj, frame)) 
AtomsBase.cell_vectors(trj::Trajectory, frame) = cell_vectors(cell(trj, frame)) 

function AtomsBase.position(trj::Trajectory, atom::Int, frame::Int)
    @argcheck 1 <= frame <= length(trj)
    @argcheck 1 <= atom <= trj.natoms
    k = (frame-1)*trj.natoms + atom
    return trj.position[k]
end

##

