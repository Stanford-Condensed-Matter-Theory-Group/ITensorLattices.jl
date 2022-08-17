
export build_graph, build_lattice, get_coords, nsites
export LatticeParams, ChainLatticeParams, TriangularLatticeParams, SquareLatticeParams, HoneycombLatticeParams

"""
Aliases
"""
nsites = nnodes
const LatticeParams = GraphParams
const ChainLatticeParams = ChainGraphParams
const TriangularLatticeParams = TriangularGraphParams
const SquareLatticeParams = SquareGraphParams
const HoneycombLatticeParams = HoneycombGraphParams

"""
Build graph from lattice. This can be useful for analyzing the bond structure
with graph algorithms and visualization of the lattice as a sanity check
"""
function build_graph(lattice::Lattice)::SimpleGraph
    edges = map(bond -> Edge(bond.s1, bond.s2), lattice)
    return SimpleGraph(edges)
end

"""
Build ITensor Lattice from undirected graph
"""
function build_lattice(graph::Graph, coords::Union{Vector{<:Coord}, Nothing}=nothing)::Lattice
    @assert !is_directed(graph) "Graph must be undirected"
    es = edges(graph)

    function buildbond(edge)
        idx1 = src(edge)
        idx2 = dst(edge)
        (x1, y1), (x2, y2) = isnothing(coords) ? ((0, 0), (0, 0)) : (coords[idx1], coords[idx2])
        return LatticeBond(idx1, idx2, x1, y1, x2, y2)
    end

    return buildbond.(es)
end

"""
Builds a dictionary mapping the site index to cartesian coordinates
"""
function get_coords(lattice::Lattice)::Vector{Coord}
    isempty(lattice) && return []

    # make index => (x, y) pairs
    getcoord(bond) = [bond.s1 => (bond.x1, bond.y1), bond.s2 => (bond.x2, bond.y2)]
    allcoords = vcat(getcoord.(lattice)...)

    # create dictionary to remove duplicates
    dict = Dict(allcoords)

    # create array or coordinates indexed by the site index
    indices = keys(dict)
    maxidx = maximum(indices)
    coordlist = Vector{Coord}(undef, maxidx)
    for i in 1:maxidx
        coordlist[i] = get(dict, i, (0, 0))
    end

    return coordlist
end

function build_lattice(params::LatticeParams; neighbor::Int=1)
    graph, coords = build_graph(params)
    neighbor_graph = get_neighbor_graph(graph, neighbor)
    return build_lattice(neighbor_graph, coords)
end
