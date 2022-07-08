
export build_graph, build_lattice, get_coords, build_chain_lattice,
    build_triangular_lattice, build_square_lattice, build_honeycomb_lattice

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

TODO: handle xy coordinates
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

############################################################################################
### Lattices
############################################################################################

function build_chain_lattice(n::Int; periodic = false, neighbor = 1)
    g, coords = build_chain_graph(n; periodic)
    graph = get_neighbor_graph(g, neighbor)
    return build_lattice(graph, coords)
end

function build_triangular_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    g, coords = build_triangular_graph(nx, ny; yperiodic)
    graph = get_neighbor_graph(g, neighbor)
    return build_lattice(graph, coords)
end

function build_square_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    g, coords = build_square_graph(nx, ny; yperiodic)
    graph = get_neighbor_graph(g, neighbor)
    return build_lattice(graph, coords)
end

"""
armchair configuration
"""
function build_honeycomb_lattice(nx::Int, ny::Int; yperiodic::Bool=false, neighbor::Int=1)::Lattice
    g, coords = build_honeycomb_graph(nx, ny; yperiodic)
    graph = get_neighbor_graph(g, neighbor)
    return build_lattice(graph, coords)
end
