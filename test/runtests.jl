using AtomsBase
using AtomsBaseTesting
using SimpleAtomsStructures
using Rotations
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
        @test all( mass(sys, :) .== mass( ref.system, :) )
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
        @test isa(sys, SimpleAtomsStructures.GeneralSystem)
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
    @testset "Utils" begin
        sys = SimpleSystem(ref.system)
        # rotation tests
        q = rand(QuatRotation)
        sys2 = sys * q
        @test all( i-> position(sys2, i) ≈ q * position(sys, i), 1:length(sys) )
        @test angle(sys, 1, 2, 3) ≈ angle(sys2, 1, 2, 3)
        @test angled(sys, 1, 2, 3) ≈ angled(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test dihedral_angled(sys, 1,2,3,4) ≈ dihedral_angled(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ q * distance_vector(sys, 1, 2)

        # translation tests
        cms = center_of_mass(sys)
        sys2 = sys - cms
        @test all( i-> position(sys2, i) ≈ position(sys, i) - cms, 1:length(sys) )
        @test angle(sys, 1, 2, 3) ≈ angle(sys2, 1, 2, 3)
        @test angled(sys, 1, 2, 3) ≈ angled(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test dihedral_angled(sys, 1,2,3,4) ≈ dihedral_angled(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ distance_vector(sys, 1, 2)

        
        sys3 = sys + sys2
        @test length(sys3) == length(sys) + length(sys2)
        
        # now with cell
        sys = GenericSystem(ref.system)

        fp = fractional_coordinates(sys, :)
        @test length(fp) == length(sys)
        fpm = fractional_coordinates_as_matrix(sys, :)
        @test size(fpm) == (3, length(sys))

        clm = cell_matrix(sys)
        clv = cell_vectors(sys)
        @test size(clm) == (3,3)
        @test all( x -> all(x[1] .≈ x[2]), zip(clv, eachcol(clm)) )
        icell = inv_cell(sys)
        tmp = icell * clm
    end
    @testset "SimpleAtom" begin
        sys = GenericSystem(ref.system)
        va = sys[:]
        @test all( k -> k in atomkeys(sys), atomkeys(va) )
        @test all( k -> k in atomkeys(va), atomkeys(sys) )
        @test all( mass(va, :) .≈ mass(sys, :) )
        @test all( species(va, :) .== species(sys, :) )
        @test all( position(va, :) .≈ position(sys, :) )
        @test all( velocity(va, :) .≈ velocity(sys, :) )
    end
end
