

struct SystemView{D, LU, TT} <: AtomsBase.AbstractSystem{D}
    parent_trajectory::TT
    frame::Int
    function SystemView(trj::AbstractTrajectory{D, LU}, frame::Int) where{D, LU}
        @argcheck 1 <= frame <= length(trj)
        new{D, LU, typeof(trj)}(trj, frame)
    end
end

function Base.getindex(sv::SystemView, i::Int)
    return SimpleAtom(species(sv, i), position(sv, i))
end

function Base.length(sv::SystemView)
    return natoms(sv.parent_trajectory)
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

AtomsBase.atomkeys(sv::SystemView) = atomkeys(sv.parent_trajectory)
AtomsBase.hasatomkey(sv::SystemView, key::Symbol) = in(key, atomkeys(sv) )
AtomsBase.cell(sv::SystemView) = cell(sv.parent_trajectory, sv.frame)

AtomsBase.species(sv::SystemView, i) = species(sv.parent_trajectory, i)
AtomsBase.mass(sv::SystemView, i) = mass(sv.parent_trajectory, i)

AtomsBase.position(sv::SystemView, i) = position(sv.parent_trajectory, i, sv.frame)
AtomsBase.velocity(sv::SystemView, i) = velocity(sv.parent_trajectory, i, sv.frame)



## struct Definitions

mutable struct Trajectory{D, LU, TP, TC} <: AbstractTrajectory{D, LU}
    constants::TC
    position::Vector{SVector{D, TP}}
    function Trajectory(
        spc::AbstractVector{ChemicalSpecies},
        pos::AbstractVector{SVector{D, TP}};
        cell=IsolatedCell(D, TP)
    ) where {D, TP<:Unitful.Length}
        @argcheck length(spc) == length(pos)
        LU = unit(TP)
        # need copies to not break things - trajectory is mutable
        tmp = ( species=deepcopy(spc), natoms=length(spc), cell=cell )
        new{D, LU, TP, typeof(tmp)}( tmp, deepcopy(pos))
    end    
end


mutable struct VelocityTrajectory{D, LU, TP, TV, TC} <: AbstractTrajectory{D, LU}
    constants::TC
    position::Vector{SVector{D, TP}}
    velocity::Vector{SVector{D, TV}}
    function VelocityTrajectory(
        spc::AbstractVector{ChemicalSpecies},
        pos::AbstractVector{SVector{D, TP}},
        vel::AbstractVector{SVector{D, TV}};
        cell::Union{IsolatedCell, PeriodicCell}=IsolatedCell(D, TP)
    ) where {D, TP, TV}
        @argcheck length(spc) == length(pos) == length(vel)
        LU = unit(TP)
        tmp = ( species=deepcopy(spc), natoms=length(spc), cell=cell )
        new{D, LU, TP, TV, typeof(tmp)}(
            tmp,
            deepcopy(pos),
            deepcopy(vel),
        )
    end    
end


## Constructors


function Trajectory(sys::AbstractSystem)
    spc = species(sys, :) 
    pos = position(sys, :) 
    return Trajectory(spc, pos; cell=cell(sys))
end

function Trajectory(vsys::AbstractVector{<:AbstractSystem})
    tmp = Trajectory(vsys[1])
    append!(tmp, vsys[2:end])
    return tmp
end

function VelocityTrajectory(sys::AbstractSystem)
    spc = species(sys, :) 
    pos = position(sys, :)
    vel = velocity(sys, :)
    return VelocityTrajectory(spc, pos, vel; cell=cell(sys))
end

function VelocityTrajectory(vsys::AbstractVector{<:AbstractSystem})
    tmp = VelocityTrajectory(vsys[1])
    append!(tmp, vsys[2:end])
    return tmp
end

##

function Base.length(trj::AbstractTrajectory)
    return length(trj.position) // trj.constants.natoms |> Int
end

@inline natoms(trj::AbstractTrajectory) = trj.constants.natoms

Base.eltype(::Type{Trajectory{D, LU, TP, TC}}) where{D, LU, TP, TC} =
    SystemView{D, LU, Trajectory{D, LU, TP, TC}}
Base.eltype(::Type{VelocityTrajectory{D, LU, TP, TV, TC}}) where{D, LU, TP, TV, TC} =
    SystemView{D, LU, VelocityTrajectory{D, LU, TP, TV, TC}}
Base.size(trj::AbstractTrajectory) = (length(trj), )

Base.show(io::IO, trj::Trajectory) =
    print(io, "Trajectory with ", length(trj), " frames of ", natoms(trj), " atoms")
Base.show(io::IO, trj::VelocityTrajectory) =
    print(io, "VelocityTrajectory with ", length(trj), " frames of ", natoms(trj), " atoms")
Base.show(io::IO, ::MIME"text/plain", trj::AbstractTrajectory) = show(io, trj)


function Base.push!(trj::Trajectory{D}, pos::AbstractVector{SVector{D, TP}}) where{D, TP<:Unitful.Length}
    @argcheck length(pos) == natoms(trj)
    append!(trj.position, pos)
    return trj
end

function Base.push!(
    trj::VelocityTrajectory{D}, 
    pos::AbstractVector{SVector{D, TP}}, 
    vel::AbstractVector{SVector{D, TV}}
) where{D, TP, TV}
    @argcheck length(pos) == length(vel) == natoms(trj)
    append!(trj.position, pos)
    append!(trj.velocity, vel)
    return trj
end

function Base.push!(trj::Trajectory{D}, sys::AbstractSystem{D}; check_species=false) where{D}
    @argcheck trj.natoms == length(sys)
    @argcheck cell(trj) == cell(sys)
    if check_species
        @argcheck all( i -> species(trj, i) === species(sys, i) for i = 1:natoms(trj) )
    end
    append!(trj.position, position(sys, :))
    return trj
end

function Base.push!(trj1::Trajectory{D}, trj2::AbstractTrajectory{D}; check_species=false) where{D}
    @argcheck natoms(trj1) == natoms(trj2)
    @argcheck cell(trj1) == cell(trj2)
    if check_species
        @argcheck all( i -> species(trj1, i) === species(trj2, i) for i = 1:natoms(trj1) )
    end
    append!(trj1.position, trj2.position)
    return trj1
end

function Base.push!(trj::VelocityTrajectory{D}, sys::AbstractSystem{D}; check_species=false) where{D}
    @argcheck trj.natoms == length(sys)
    @argcheck cell(trj) == cell(sys)
    if check_species
        @argcheck all( i -> species(trj, i) === species(sys, i) for i = 1:natoms(trj) )
    end
    append!(trj.position, position(sys, :))
    append!(trj.velocity, velocity(sys, :))
    return trj  
end

function Base.push!(trj1::VelocityTrajectory{D}, trj2::VelocityTrajectory{D}; check_species=false) where{D}
    @argcheck natoms(trj1) == natoms(trj2)
    @argcheck cell(trj1) == cell(trj2)
    if check_species
        @argcheck all( i -> species(trj1, i) === species(trj2, i) for i = 1:natoms(trj1) )
    end
    append!(trj1.position, trj2.position)
    append!(trj1.velocity, trj2.velocity)
    return trj1
end

function Base.append!(trj::Trajectory{D}, addtraj::AbstractVector{<:AbstractSystem{D}}) where{D}
    # NOTE species are not tested
    @argcheck all( natoms(trj) == length(sys) for sys in addtraj )
    pos = mapreduce(vcat, addtraj) do frame
        position(frame, :)
    end
    append!(trj.position, pos)
    return trj
end

function Base.append!(trj::VelocityTrajectory{D}, addtraj::AbstractVector{<:AbstractSystem{D}}) where{D}
    # NOTE species are not tested
    @argcheck all( natoms(trj) == length(sys) for sys in addtraj )
    pos = mapreduce(vcat, addtraj) do frame
        position(frame, :)
    end
    vel = mapreduce(vcat, addtraj) do frame
        velocity(frame, :)
    end
    append!(trj.position, pos)
    append!(trj.velocity, vel)
    return trj
end

function Base.append!(trj1::Trajectory{D}, trj2::AbstractTrajectory{D}) where{D}
    return push!(trj1, trj2)
end

function Base.append!(trj1::VelocityTrajectory{D}, trj2::VelocityTrajectory{D}) where{D}
    return push!(trj1, trj2)
end

function Base.getindex(trj::AbstractTrajectory, i::Int)
    @argcheck 1 <= i <= length(trj)
    return SystemView(trj, i)
end

Base.haskey(trj::Trajectory, x::Symbol) = in(x, keys(trj) )
Base.keys(::Trajectory) = (:cell_vectors, :periodicity)

AtomsBase.atomkeys(::Trajectory) = (:position, :species)
AtomsBase.hasatomkey(trj::AbstractTrajectory, key) = key in atomkeys(trj)
AtomsBase.atomkeys(::VelocityTrajectory) = (:position, :velocity, :species)


AtomsBase.species(trj::AbstractTrajectory, i, frame) = species(trj, i)
AtomsBase.species(trj::AbstractTrajectory, i) = view(trj.constants.species, i)
AtomsBase.species(trj::AbstractTrajectory, ::Colon) = trj.constants.species
AtomsBase.species(tra::AbstractTrajectory, i::Int) = tra.constants.species[i]

AtomsBase.mass(trj::AbstractTrajectory, i, frame) = mass(trj, i)
AtomsBase.mass(trj::AbstractTrajectory, i) = mass.( species(trj, i) )
AtomsBase.cell(trj::AbstractTrajectory) = trj.constants.cell
AtomsBase.cell(trj::AbstractTrajectory, ::Int) = trj.constants.cell
AtomsBase.periodicity(trj::AbstractTrajectory, frame) = periodicity(cell(trj, frame)) 
AtomsBase.cell_vectors(trj::AbstractTrajectory, frame) = cell_vectors(cell(trj, frame)) 

function AtomsBase.position(trj::AbstractTrajectory, atom::Int, frame::Int)
    @argcheck 1 <= frame <= length(trj)
    @argcheck 1 <= atom <= natoms(trj)
    k = (frame-1)*natoms(trj) + atom
    return trj.position[k]
end

function AtomsBase.position(trj::AbstractTrajectory, ::Colon, frame::Int)
    @argcheck 1 <= frame <= length(trj)
    k = (frame-1)*natoms(trj) + 1 : frame*natoms(trj)
    return view(trj.position, k)
end

function AtomsBase.position(trj::AbstractTrajectory, atoms, frame)
    pos = reshape(trj.position, natoms(trj), length(trj))
    return view(pos, atoms, frame)
end

function AtomsBase.velocity(trj::VelocityTrajectory, atom::Int, frame::Int)
    @argcheck 1 <= frame <= length(trj)
    @argcheck 1 <= atom <= natoms(trj)
    k = (frame-1)*natoms(trj) + atom
    return trj.velocity[k]
end

function AtomsBase.velocity(trj::VelocityTrajectory, ::Colon, frame::Int)
    @argcheck 1 <= frame <= length(trj)
    k = (frame-1)*natoms(trj) + 1 : frame*natoms(trj)
    return view(trj.velocity, k)
end

function AtomsBase.velocity(trj::VelocityTrajectory, atoms, frame)
    vel = reshape(trj.velocity, natoms(trj), length(trj))
    return view(vel, atoms, frame)
end

##
