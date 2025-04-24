module AtomsTrajectories

using ArgCheck
using AtomsBase
using LinearAlgebra: norm
using StaticArrays
using Unitful
import ..SimpleAtomsStructures

export SystemView, SimpleTrajectory, VelocityTrajectory

abstract type AbstractTrajectory{D, LU} <: AbstractVector{AbstractSystem{D}} end

include("trajectory.jl")
include("general_methods.jl")

end
