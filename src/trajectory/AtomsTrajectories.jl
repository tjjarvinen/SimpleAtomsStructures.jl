module AtomsTrajectories

using ArgCheck
using AtomsBase
using Unitful
using StaticArrays

export SystemView, SimpleTrajectory, VelocityTrajectory

abstract type AbstractTrajectory{D, LU} <: AbstractVector{AbstractSystem{D}} end

include("trajectory.jl")

end
