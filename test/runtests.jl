using AtomsBase
using AtomsBaseTesting
using SimpleAtomsStructures
using Test

@testset "SimpleAtomsStructures.jl" begin
    # Write your tests here.
    ref = make_test_system()
    @testset "SimpleSystem" begin
        sys = SimpleSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy] 
    end
    @testset "SimpleVelocitySystem" begin
        sys = SimpleVelocitySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
    end
    @testset "AtomicPropertySystem" begin
        sys = AtomicPropertySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "CellSystemSystem" begin
        sys = CellSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( sys[:periodicity] .== ref.periodicity )
        @test all( sys[:cell_vectors] .≈ ref.cell_vectors )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "GeneralSystem" begin
        sys = GenericSystem(ref.system)
        @test isa(sys, GeneralSystem)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
        @test all( x-> haskey(sys, x), keys(ref.sysprop) )
        @test all( x-> haskey(ref.sysprop, x), keys(sys) )
        @test all( [all(sys[k] .== v) for (k,v) in pairs(ref.sysprop)] )
        @test_throws KeyError sys[:dummy]
    end
end
