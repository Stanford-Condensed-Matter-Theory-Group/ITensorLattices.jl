using ITensorLattices
using Test

@testset "Lattice Tests" begin
    params = ChainLatticeParams(N = 5, periodic = true)
    n = nsites(params)
    @test n == 5
    l = build_lattice(params)
    @test length(l) == 5

    params = SquareLatticeParams(Nx = 10, Ny = 3)
    n = nsites(params)
    @test n == 30
    l = build_lattice(params)
    @test length(l) == 47

    params = TriangularLatticeParams(Nx = 20, Ny = 30)
    n = nsites(params)
    @test n == 600
    l = build_lattice(params)
    @test length(l) == 1701

    params = HoneycombLatticeParams(Nx = 100, Ny = 100, yperiodic = true)
    n = nsites(params)
    @test n == 10000
    l = build_lattice(params)
    @test length(l) == 14900

    # neighbors
    params = HoneycombLatticeParams(Nx = 20, Ny = 4, yperiodic = true)
    l = build_lattice(params, neighbor=2)
    @test length(l) == 224
end
