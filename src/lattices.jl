
export build_chain_lattice, build_triangular_lattice, build_square_lattice, build_honeycomb_lattice

function build_chain_lattice(n::Int; periodic = false, neighbor = 1)
    nbonds = n - (periodic ? 0 : neighbor)
    chain = Lattice(undef, nbonds)

    for i in 1:nbonds
        nexti = ((i + neighbor - 1) % n) + 1
        chain[i] = LatticeBond(i, nexti, i, 0, nexti, 0)
    end

    return chain
end

function build_triangular_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    supported_neighbors = [1, 2]

    if neighbor == 1
        # use ITensors lattice
        return triangular_lattice(nx, ny; yperiodic)
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
        return square_lattice(nx, ny; yperiodic)
    elseif neighbor == 2
        error("Not implemented!")
    else
        error("only supports $supported_neighbors neighbors")
    end
end

"""
Build honeycomb lattice in armchair configuration

y = 3 o - o - o - o - o - o - o - o
        /      \\       /       \
y = 2 o - o - o - o - o - o - o - o
        /      \\       /       \
y = 1 o - o - o - o - o - o - o - o
  x = 1   2   3   4   5   6   7   8
"""
function build_honeycomb_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    # only handle nearest neighbor
    @assert neighbor == 1 "Only handles neighbor = 1"

    getidx(x, y) = (x - 1) * ny + y
    wrap(x, k) = mod1(mod1(x, k) + k, k)

    function getcoords(x, y)
        xcoord = (x - 1) * 3/4 + ((x % 2) == 0 ? -1/4 : 0)
        ycoord = (y - 1) * sqrt(3) + ((x ÷ 2) % 2 == 0 ? sqrt(3) / 2 : 0)
        return xcoord, ycoord
    end

    lattice = Lattice(undef, 0)

    function bond!(x1, y1, x2, y2)
        idx1 = getidx(x1, y1)
        idx2 = getidx(x2, y2)
        x1coord, y1coord = getcoords(x1, y1)
        x2coord, y2coord = getcoords(x2, y2)
        push!(lattice, LatticeBond(idx1, idx2, x1coord, y1coord, x2coord, y2coord))
    end

    for x in 1:(nx - 1)
        for y in 1:ny
            # bind to next row
            bond!(x, y, x + 1, y)

            # Alternate between linking up and linking down on odd rows
            if x % 2 == 1
                crossy = y + ((x ÷ 2) % 2 == 0 ? 1 : -1)
                if crossy ∈ 1:ny
                    bond!(x, y, x + 1, crossy)
                elseif yperiodic
                    bond!(x, y, x + 1, wrap(crossy, ny))
                end
            end
        end
    end

    return lattice
end
