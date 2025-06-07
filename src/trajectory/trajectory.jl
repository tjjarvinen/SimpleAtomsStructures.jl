struct SystemView{D, LU, TT} <: AtomsBase.AbstractSystem{D}
    parent_trajectory::TT
    frame::Int
    function SystemView(trj::AbstractTrajectory{D, LU}, frame::Int) where{D, LU}
        @argcheck 1 <= frame <= length(trj)
        new{D, LU, typeof(trj)}(trj, frame)
    end
end

function Base.getindex(sv::SystemView, i::Int)
    return _get_atom(sv.parent_trajectory, i, sv.frame)
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

mutable struct SimpleTrajectory{D, LU, TP, TC} <: AbstractTrajectory{D, LU}
    constants::TC
    position::Vector{SVector{D, TP}}
    function SimpleTrajectory(
        spc::AbstractVector{ChemicalSpecies},
        pos::AbstractVector{SVector{D, TP}};
        cell=IsolatedCell(D, TP)
    ) where {D, TP<:Unitful.Length}
        @argcheck length(spc) == length(pos)
        LU = unit(TP)
        # need copies to not break things - trajectory is mutable
        tmp = ( atom_constants= (species=deepcopy(spc), ), cell=cell, natoms=length(spc))
        new{D, LU, TP, typeof(tmp)}( tmp, deepcopy(pos) )
    end    
end


mutable struct VelocityTrajectory{D, LU, TP, TC, TV} <: AbstractTrajectory{D, LU}
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
        tmp = ( atom_constants= (species=deepcopy(spc), ), cell=cell, natoms=length(spc))
        new{D, LU, TP, typeof(tmp), TV}(
            tmp,
            deepcopy(pos),
            deepcopy(vel),
        )
    end    
end


## Constructors


function SimpleTrajectory(sys::AbstractSystem)
    spc = species(sys, :) 
    pos = position(sys, :) 
    return SimpleTrajectory(spc, pos; cell=cell(sys))
end

function SimpleTrajectory(vsys::AbstractVector{<:AbstractSystem})
    tmp = SimpleTrajectory(vsys[1])
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

function _get_atom(traj::SimpleTrajectory, i::Int, frame::Int)
    return SimpleAtomsStructures.SimpleAtom(species(traj, i, frame), position(traj, i, frame))
end

function _get_atom(traj::VelocityTrajectory, i::Int, frame::Int)
    return SimpleAtomsStructures.SimpleAtom(
        species(traj, i, frame), 
        position(traj, i, frame), 
        velocity(traj, i, frame)
    )
end


function Base.length(trj::AbstractTrajectory)
    return length(trj.position) // trj.constants.natoms |> Int
end

@inline natoms(trj::AbstractTrajectory) = trj.constants.natoms

Base.eltype(::Type{SimpleTrajectory{D, LU, TP, TC}}) where{D, LU, TP, TC} =
    SystemView{D, LU, SimpleTrajectory{D, LU, TP, TC}}
Base.eltype(::Type{VelocityTrajectory{D, LU, TP, TV, TC}}) where{D, LU, TP, TV, TC} =
    SystemView{D, LU, VelocityTrajectory{D, LU, TP, TV, TC}}
Base.size(trj::AbstractTrajectory) = (length(trj), )

Base.show(io::IO, trj::SimpleTrajectory) =
    print(io, "SimpleTrajectory with ", length(trj), " frames of ", natoms(trj), " atoms")
Base.show(io::IO, trj::VelocityTrajectory) =
    print(io, "VelocityTrajectory with ", length(trj), " frames of ", natoms(trj), " atoms")
Base.show(io::IO, ::MIME"text/plain", trj::AbstractTrajectory) = show(io, trj)


function Base.push!(trj::SimpleTrajectory{D}, pos::AbstractVector{SVector{D, TP}}) where{D, TP<:Unitful.Length}
    @argcheck length(pos) == natoms(trj)
    append!(trj.position, pos)
    return trj
end

function Base.push!(trj::SimpleTrajectory{D}, pos::AbstractMatrix{TP}) where{D, TP<:Unitful.Length}
    @argcheck size(pos, 1) == D
    @argcheck size(pos, 2) == natoms(trj)
    tmp = reinterpret(reshape, SVector{D, TP}, pos)
    push!(trj, tmp)
    return trj
end

function Base.push!(
    trj::VelocityTrajectory{D}, 
    pos::AbstractVector{SVector{D, TP}}, 
    vel::AbstractVector{SVector{D, TV}}
) where{D, TP<:Unitful.Length, TV<:Unitful.Velocity}
    @argcheck length(pos) == length(vel) == natoms(trj)
    append!(trj.position, pos)
    append!(trj.velocity, vel)
    return trj
end

function Base.push!(
    trj::VelocityTrajectory{D}, 
    pos::AbstractMatrix{TP}, 
    vel::AbstractMatrix{TV}
) where{D, TP<:Unitful.Length, TV<:Unitful.Velocity}
    @argcheck size(pos) == size(vel)
    @argcheck size(pos, 1) == D
    @argcheck size(pos, 2) == natoms(trj)
    tmp = reinterpret(reshape, SVector{D, TP}, pos)
    tmp2 = reinterpret(reshape, SVector{D, TV}, vel)
    push!(trj, tmp, tmp2)
    return trj
end

function Base.push!(trj::SimpleTrajectory{D}, sys::AbstractSystem{D}; check_species=false) where{D}
    @argcheck trj.natoms == length(sys)
    @argcheck cell(trj) == cell(sys)
    if check_species
        @argcheck all( i -> species(trj, i) === species(sys, i) for i = 1:natoms(trj) )
    end
    append!(trj.position, position(sys, :))
    return trj
end

function Base.push!(trj1::SimpleTrajectory{D}, trj2::AbstractTrajectory{D}; check_species=false) where{D}
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

function Base.append!(trj::SimpleTrajectory{D}, addtraj::AbstractVector{<:AbstractSystem{D}}) where{D}
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

function Base.append!(trj1::SimpleTrajectory{D}, trj2::AbstractTrajectory{D}) where{D}
    return push!(trj1, trj2)
end

function Base.append!(trj1::VelocityTrajectory{D}, trj2::VelocityTrajectory{D}) where{D}
    return push!(trj1, trj2)
end

function Base.getindex(trj::AbstractTrajectory, i::Int)
    @argcheck 1 <= i <= length(trj)
    return SystemView(trj, i)
end

Base.haskey(trj::SimpleTrajectory, x::Symbol) = in(x, keys(trj) )
Base.keys(::SimpleTrajectory) = (:cell_vectors, :periodicity)

AtomsBase.atomkeys(trj::SimpleTrajectory) = (:position, :species, keys(trj.constants.atom_constants)...)
AtomsBase.hasatomkey(trj::AbstractTrajectory, key) = key in atomkeys(trj)
AtomsBase.atomkeys(trj::VelocityTrajectory) = (:position, :velocity, :species, keys(trj.constants.atom_constants)...)


AtomsBase.species(trj::AbstractTrajectory, i, frame) = species(trj, i)
AtomsBase.species(trj::AbstractTrajectory, i) = view(trj.constants.atom_constants.species, i)
AtomsBase.species(trj::AbstractTrajectory, ::Colon) = trj.constants.atom_constants.species
AtomsBase.species(tra::AbstractTrajectory, i::Int) = tra.constants.atom_constants.species[i]

AtomsBase.mass(trj::AbstractTrajectory, i, frame) = mass(trj, i)
function AtomsBase.mass(trj::AbstractTrajectory, i)
    if haskey(trj.constants.atom_constants, :mass)
        @view trj.constants.atom_constants.mass[i]
    else
        mass.( species(trj, i) )
    end
end

function AtomsBase.mass(trj::AbstractTrajectory, i::Int)
    if haskey(trj.constants.atom_constants, :mass)
        trj.constants.atom_constants.mass[i]
    else
        mass( species(trj, i) )
    end
end

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

function AtomsBase.set_species!(trj::AbstractTrajectory, i::Int, spc::ChemicalSpecies)
    @argcheck 1 <= i <= natoms(trj)
    trj.constants.species[i] = spc
    return trj 
end

function AtomsBase.set_species!(trj::AbstractTrajectory, i, spc::AbstractVector{ChemicalSpecies})
    trj.constants.species[i] .= spc
    return trj
end

function AtomsBase.set_species!(trj::AbstractTrajectory, spc::AbstractVector{ChemicalSpecies})
    @argcheck length(spc) == natoms(trj)
    trj.constants.species .= spc
    return trj
end

##


function SimpleAtomsStructures.fractional_coordinates(traj::AbstractTrajectory, atom, frame)
    pos = position(traj, atom, frame)
    ce = cell(traj)
    tmp = reinterpret(reshape, (eltype∘eltype)(pos), vec(pos))
    tmp = SimpleAtomsStructures.inv_cell(ce) * tmp
    tmp = reinterpret(reshape, SVector{3, eltype(tmp)}, tmp)
    return reshape(tmp, size(pos))
end

function SimpleAtomsStructures.fractional_coordinates(traj::AbstractTrajectory, atom::Int, frame::Int)
    SimpleAtomsStructures.fractional_coordinates(traj[frame], atom)
end

function _distance_vector(traj::AbstractTrajectory, atom1::Int, atom2::Int, frame)
    r = position(traj, atom1, frame)
    r2 = position(traj, atom2, frame)
    tmp = r2 .- r 
    tmp = SimpleAtomsStructures.distance_vector(cell(traj), vec(tmp))
    return reshape(tmp, size(r)) 
end

function SimpleAtomsStructures.distance_vector(traj::AbstractTrajectory, atom1::Int, atom2::Int, ::Int)
    _distance_vector(traj, atom1, atom2, 1)
end

function SimpleAtomsStructures.distance_vector(traj::AbstractTrajectory, atom1::Int, atom2::Int, frame)
    _distance_vector(traj, atom1, atom2, frame)
end

function SimpleAtomsStructures.distance_vector(traj::AbstractTrajectory, atom1::Int, atom2::Int, ::Colon)
    r1 = position(traj, atom1, :)
    r2 = position(traj, atom2, :)
    Δr = r2 .- r1
    return SimpleAtomsStructures.distance_vector(cell(traj), vec(Δr))
end


function _distance(traj::AbstractTrajectory, atom1, atom2, frame)
    tmp = SimpleAtomsStructures.distance_vector(traj, atom1, atom2, frame)
    return _extract_norm(tmp)
end

function SimpleAtomsStructures.distance(traj::AbstractTrajectory, atom1::Int, atom2::Int, frame::Int)
    _distance(traj, atom1, atom2, frame)
end

function SimpleAtomsStructures.distance(traj::AbstractTrajectory, atom1::Int, atom2::Int, frame::Colon)
    _distance(traj, atom1, atom2, frame)
end

function SimpleAtomsStructures.distance(traj::AbstractTrajectory, atom1::Int, atom2::Int, frames)
    _distance(traj, atom1, atom2, frames)
end


_extract_norm(x::SVector) = norm(x)
_extract_norm(x::AbstractVector{<:SVector}) = norm.(x)