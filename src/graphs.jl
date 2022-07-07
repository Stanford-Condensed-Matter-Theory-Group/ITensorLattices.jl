
export build_graph, build_lattice

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
function build_lattice(graph::Graph)::Lattice
    @assert !is_directed(graph) "Graph must be undirected"
    es = edges(graph)
    latt = map(edge -> LatticeBond(src(edge), dst(edge)), es)
end

"""
Builds a dictionary mapping the site index to cartesian coordinates
"""
function get_coords(lattice::Lattice)
    # make index => (x, y) pairs
    getcoords(bond) = [bond.s1 => (bond.x1, bond.y1), bond.s2 => (bond.x2, bond.y2)]
    allcoords = vcat(getcoords.(lattice)...)

    # create dictionary to remove duplicates
    dict = Dict(allcoords)

    # create array or coordinates indexed by the site index
    indices = keys(dict)
    maxidx = maximum(indices)
    coordlist = Vector(undef, maxidx)
    for (idx, coords) in dict
        coordlist[idx] = coords
    end

    return coordlist
end
