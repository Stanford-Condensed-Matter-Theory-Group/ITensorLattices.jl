
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
