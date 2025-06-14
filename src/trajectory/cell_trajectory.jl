
mutable struct VariableVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}
    base_trajectory::TB
    cell::Vector{PeriodicCell{D, TP}}
    function VariableVolumeTrajectory(
        traj::AbstractSimpleTrajectory{D, LU, TP},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP}
        new{D, LU, TP, typeof(traj)}(
            traj,
            [cell]
        )
    end
end

function VariableVolumeTrajectory(sys)
    base = SimpleVelocityTrajectory(sys)
    if ! (cell(sys) isa PeriodicCell)
        return base
    end
    return VariableVolumeTrajectory(base, cell(sys))
end

@inline n_atoms(sys::AbstractCellTrajectory) = n_atoms(sys.base_trajectory)

Base.size(sys::AbstractCellTrajectory) = size(sys.base_trajectory)

function Base.getindex(sys::VariableVolumeTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(sys) BoundsError(sys, i)
    return SimpleAtomsStructures.CellSystemView(
        sys.base_trajectory[i],
        sys.cell[i]
    )
end

function Base.append!(traj::VariableVolumeTrajectory{D, LU, TP}, sys::AbstractSystem{D}) where {D, LU, TP}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    if cell(sys) isa PeriodicCell{D, TP}
        append!(traj.base_trajectory, sys)
        push!(traj.cell, cell(sys))
    else
        error("System must have a periodic cell to append to a VariableVolumeTrajectory") 
    end
    return traj
end


function Base.eltype(::AbstractCellTrajectory{D, LU}) where {D, LU}
    return SimpleAtomsStructures.CellSystemView{D, LU}
end


##

mutable struct ConstantVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}
    base_trajectory::TB
    cell::PeriodicCell{D, TP}
    function ConstantVolumeTrajectory(
        traj::AbstractSimpleTrajectory{D, LU, TP},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP}
        new{D, LU, TP, typeof(traj)}(
            traj,
            cell
        )
    end
end

function ConstantVolumeTrajectory(sys)
    base = SimpleVelocityTrajectory(sys)
    if ! (cell(sys) isa PeriodicCell)
        return base
    end
    return ConstantVolumeTrajectory(base, cell(sys))
end

function Base.getindex(sys::ConstantVolumeTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(sys) BoundsError(sys, i)
    return SimpleAtomsStructures.CellSystemView(
        sys.base_trajectory[i],
        sys.cell
    )   
end


function Base.append!(traj::ConstantVolumeTrajectory{D, LU, TP}, sys::AbstractSystem{D}) where {D, LU, TP}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    if cell(sys) == traj.cell
        append!(traj.base_trajectory, sys)
    else
        error("System cell must match the trajectory cell to append") 
    end
    return traj
end