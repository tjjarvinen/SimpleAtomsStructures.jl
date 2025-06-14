
mutable struct SimpleTrajectory{D, LU, TP} <: AbstractSimpleTrajectory{D, LU, TP}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    n_atoms::Int
    function SimpleTrajectory(
        species::AbstractVector{ChemicalSpecies},
        positions::AbstractVector{SVector{D, TP}}
    ) where {D, TP<:Unitful.Length}
        @argcheck length(species) > 0 "Species vector must not be empty"
        n_atoms = length(species)
        @argcheck length(positions) == n_atoms "Number of positions must match number of species"
        LU = unit(TP)
        new{D, LU, TP}(deepcopy(species), deepcopy(positions), n_atoms)
    end
end


function SimpleTrajectory(sys)
    return SimpleTrajectory(
        species(sys, :),
        position(sys, :)
    )
end


function Base.append!(traj::SimpleTrajectory{D}, pos::AbstractVector{SVector{D, TP}}) where {D, TP<:Unitful.Length}
    @argcheck length(pos) == n_atoms(traj) "Position vector must have the same number of atoms as the trajectory"
    append!(traj.position, pos)
    return traj
end

function Base.append!(traj::SimpleTrajectory{D}, sys::AbstractSystem{D}) where {D}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    @argcheck all(species(sys, :) .=== traj.species) "Species must match the trajectory species"
    return Base.append!(traj, position(sys, :))    
end

function Base.push!(traj::SimpleTrajectory{D}, pos::AbstractVector{SVector{D, TP}}) where {D, TP<:Unitful.Length}
    @argcheck length(pos) == n_atoms(traj) "Position vector must have the same number of atoms as the trajectory"
    append!(traj.position, pos)
    return traj
end

function Base.push!(trj::SimpleTrajectory{D}, pos::AbstractMatrix{TP}) where{D, TP<:Unitful.Length}
    @argcheck size(pos, 1) == D
    @argcheck size(pos, 2) == n_atoms(trj)
    tmp = reinterpret(reshape, SVector{D, TP}, pos)
    push!(trj, tmp)
    return trj
end


function Base.eltype(::SimpleTrajectory{D, LU, TP}) where {D, LU, TP}
    return SimpleAtomsStructures.SimpleSystemView{D, LU, TP}
end


function Base.getindex(traj::SimpleTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(traj) BoundsError(traj, i)
    k = traj.n_atoms * (i-1) +1  : traj.n_atoms * i
    return SimpleAtomsStructures.SimpleSystemView(view(traj.species, 1:traj.n_atoms), view(traj.position, k))
end

Base.size(traj::AbstractSimpleTrajectory) = (Int(length(traj.position)//traj.n_atoms), )


##

mutable struct SimpleVelocityTrajectory{D, LU, TP, TV} <: AbstractSimpleTrajectory{D, LU, TP}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    velocity::Vector{SVector{D, TV}}
    n_atoms::Int
    function SimpleVelocityTrajectory(
        species::AbstractVector{ChemicalSpecies}, 
        positions::AbstractVector{SVector{D, TP}}, 
        velocities::AbstractVector{SVector{D, TV}}
        ) where {D, TP<:Unitful.Length, TV<:Unitful.Velocity}
        @argcheck length(species) > 0 "Species vector must not be empty"
        n_atoms = length(species)
        @argcheck length(positions) == n_atoms "Number of positions must match number of species"
        @argcheck length(velocities) == n_atoms "Number of velocities must match number of species"
        LU = unit(TP)
        new{D, LU, TP, TV}(deepcopy(species), deepcopy(positions), deepcopy(velocities), n_atoms)
    end
end


function SimpleVelocityTrajectory(sys)
    if hasatomkey(sys, :velocity)
        SimpleTrajectory(sys)
    end
    return SimpleVelocityTrajectory(
        species(sys, :),
        position(sys, :),
        velocity(sys, :)
    )
end

function Base.append!(traj::SimpleVelocityTrajectory{D}, pos::AbstractVector{SVector{D, TP}}, vel::AbstractVector{SVector{D, TV}}) where {D, TP<:Unitful.Length, TV<:Unitful.Velocity}
    @argcheck length(pos) == n_atoms(traj) "Position vector must have the same number of atoms as the trajectory"
    @argcheck length(vel) == n_atoms(traj) "Velocity vector must have the same number of atoms as the trajectory"
    append!(traj.position, pos)
    append!(traj.velocity, vel)
    return traj
end

function Base.append!(traj::SimpleVelocityTrajectory{D}, sys::AbstractSystem{D}) where {D}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    @argcheck all(species(sys, :) .=== traj.species) "Species must match the trajectory species"
    return Base.append!(traj, position(sys, :), velocity(sys, :))    
end

function Base.eltype(::SimpleVelocityTrajectory{D, LU, TP, TV}) where {D, LU, TP, TV}
    return SimpleAtomsStructures.SimpleVelocitySystemView{D, LU, TP, TV}
end

function Base.getindex(traj::SimpleVelocityTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(traj) BoundsError(traj, i)
    k = traj.n_atoms * (i-1) +1  : traj.n_atoms * i
    return SimpleAtomsStructures.SimpleVelocitySystemView(
        view(traj.species, 1:traj.n_atoms), 
        view(traj.position, k), 
        view(traj.velocity, k)
    )
end
    
        
        