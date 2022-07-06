
export build_chain_lattice, build_triangular_lattice, build_square_lattice, build_honeycomb_lattice

function build_chain_lattice(n::Int; periodic = false, neighbor = 1)
    nbonds = n - (periodic ? 0 : neighbor)
    chain = Lattice(undef, nbonds)

    for i in 1:nbonds
    chain[i] = LatticeBond(i, ((i + neighbor - 1) % n) + 1)
    end

    return chain
end

function build_triangular_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    supported_neighbors = [1, 2]

    if neighbor == 1
        # use ITensors lattice
        return triangular_lattice(nx, ny, yperiodic=yperiodic)
    elseif neighbor == 2
        error("Not implemented!")
    else
        error("only supports $supported_neighbors neighbors")
    end
end

function build_square_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    supported_neighbors = [1, 2]

    if neighbor == 1
        # use ITensors lattice
        return square_lattice(nx, ny, yperiodic=yperiodic)
    elseif neighbor == 2
        error("Not implemented!")
    else
        error("only supports $supported_neighbors neighbors")
    end
end

function build_honeycomb_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    error("Not implemented!")
end
