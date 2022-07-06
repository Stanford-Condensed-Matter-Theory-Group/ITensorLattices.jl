
export visualize

"""
Visualize ITensor Lattice

TODO: explore better visualizations
"""
function visualize(lattice::Lattice)
    graph = build_graph(lattice)
    gplot(graph)
end
