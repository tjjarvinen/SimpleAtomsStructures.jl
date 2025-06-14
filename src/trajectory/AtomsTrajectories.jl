module AtomsTrajectories

using ArgCheck
using AtomsBase
using LinearAlgebra: norm
using StaticArrays
using Unitful
import ..SimpleAtomsStructures

export SimpleTrajectory
export SimpleVelocityTrajectory
export ConstantVolumeTrajectory
export VariableVolumeTrajectory
export AbstractTrajectory
export AbstractSimpleTrajectory

abstract type AbstractTrajectory{D, LU} <: AbstractVector{AbstractSystem{D}} end
abstract type AbstractSimpleTrajectory{D, LU, TP} <: AbstractTrajectory{D, LU} end
abstract type AbstractCellTrajectory{D, LU, TP} <: AbstractTrajectory{D, LU} end


@inline n_atoms(trj::AbstractTrajectory) = trj.n_atoms


include("simple_trajectory.jl")
include("cell_trajectory.jl")
include("general_methods.jl")

end
