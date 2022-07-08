
export visualize

"""
Visualize ITensor lattice
Keyword arguments are GraphRecipes arguments: https://docs.juliaplots.org/stable/generated/graph_attributes/
"""
function visualize(lattice::Lattice; use_lattice_coords=true, kwargs...)
    if isempty(lattice)
        @warn "Empty lattice — nothing to visualize"
        return
    end

    defaultargs = (
        curves=false,
        curvature_scalar=0.2,
        nodesize=0.3,
        nodeshape=:circle
    )

    graph = build_graph(lattice)

    if use_lattice_coords
        coords = get_coords(lattice)
        x = [coord[1] for coord in coords]
        y = [coord[2] for coord in coords]
        graphplot(graph; x, y, defaultargs..., kwargs...)
    else
        graphplot(graph; defaultargs..., kwargs...)
    end
end
