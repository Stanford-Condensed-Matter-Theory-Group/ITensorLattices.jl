
export build_chain_graph, build_triangular_graph, build_square_graph, build_honeycomb_graph

Coord = Tuple{<:Real, <:Real}

function get_neighbor_graph(graph, n)
    n == 1 && return graph

    A = adjacency_matrix(graph)

    # walk n steps
    N = (A ^ n) .!= 0

    # remove all neighbors that are closer than n steps away
    closer = reduce(.|, [A ^ n .!= 0 for n in 1:(n - 1)])
    N .&= .!closer

    # remove selves
    Ident = BitArray(Diagonal(ones(size(A, 1))))
    N .&= .!Ident

    return SimpleGraph(N)
end


############################################################################################
### Graph Types
############################################################################################s

function build_chain_graph(n::Int; periodic = false)
    coords = [(i - 1, 0) for i in 1:n]

    edges = Edge.([(i, i+1) for i in 1:(n-1)])
    periodic && push!(edges, Edge(n, 1))

    graph = SimpleGraph(edges)

    return graph, coords
end

function build_triangular_graph(nx::Int, ny::Int; yperiodic::Bool=false)
    yperiodic = yperiodic && (ny > 2)

    n = nx * ny
    nbonds = 3n - 2ny + (yperiodic ? 0 : -2nx + 1)
    edges = Vector{Tuple{Int, Int}}(undef, nbonds)
    coords = Vector{Coord}(undef, n)

    b = 0
    for i in 1:n
        x = (i - 1) ÷ ny + 1
        y = mod1(i, ny)

        # set coordinates
        coords[i] = (x - 1, y - 1)

        # x-direction bonds
        if x < nx
            edges[b += 1] = (i, i + ny)
        end

        # 2d bonds
        if ny > 1
            # vertical / y-periodic diagonal bond
            if (i + 1 <= n) && ((y < ny) || yperiodic)
                edges[b += 1] = (i, i + 1)
            end
            # periodic vertical bond
            if yperiodic && y == 1
                edges[b += 1] = (i, i + ny - 1)
            end
            # diagonal bonds
            if x < nx && y < ny
                edges[b += 1] = (i, i + ny + 1)
            end
        end
    end

    graph = SimpleGraph(Edge.(edges))

    return graph, coords
end

function build_square_graph(nx::Int, ny::Int; yperiodic::Bool=false)
    yperiodic = yperiodic && (ny > 2)

    n = nx * ny
    nbonds = 2n - ny + (yperiodic ? 0 : -nx)
    edges = Vector{Tuple{Int, Int}}(undef, nbonds)
    coords = Vector{Coord}(undef, n)

    b = 0
    for i in 1:n
        x = (i - 1) ÷ ny + 1
        y = mod1(i, ny)

        # set coordinates
        coords[i] = (x - 1, y - 1)

        if x < nx
            edges[b += 1] = (i, i + ny)
        end

        if ny > 1
            if y < ny
                edges[b += 1] = (i, i + 1)
            elseif yperiodic
                edges[b += 1] = (i, i - ny + 1)
            end
        end
    end

    graph = SimpleGraph(Edge.(edges))

    return graph, coords
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
function build_honeycomb_graph(nx::Int, ny::Int; yperiodic::Bool=false)
    getidx(x, y) = (x - 1) * ny + y
    wrap(x, k) = mod1(mod1(x, k) + k, k)

    # Build graph
    yperiodic = yperiodic && (ny > 2)

    edges = Vector{Tuple{Int, Int}}()

    function add_edge!(x1, y1, x2, y2)
        idx1 = getidx(x1, y1)
        idx2 = getidx(x2, y2)
        push!(edges, (idx1, idx2))
    end

    for x in 1:(nx - 1)
        for y in 1:ny
            # bind to next row
            add_edge!(x, y, x + 1, y)

            # Alternate between linking up and linking down on odd rows
            if x % 2 == 1
                crossy = y + ((x ÷ 2) % 2 == 0 ? 1 : -1)
                if crossy ∈ 1:ny
                    add_edge!(x, y, x + 1, crossy)
                elseif yperiodic
                    add_edge!(x, y, x + 1, wrap(crossy, ny))
                end
            end
        end
    end

    graph = SimpleGraph(Edge.(edges))

    # Calculate node coordinates
    function getcoords(x, y)
        xcoord = (x - 1) * 3/4 + ((x % 2) == 0 ? -1/4 : 0)
        ycoord = (y - 1) * sqrt(3) + ((x ÷ 2) % 2 == 0 ? sqrt(3)/2 : 0)
        return xcoord, ycoord
    end

    coords = Vector{Coord}(undef, nx * ny)
    for x in 1:nx
        for y in 1:ny
            xcoord = (x - 1) * 3/4 + ((x % 2) == 0 ? -1/4 : 0)
            ycoord = (y - 1) * sqrt(3) + ((x ÷ 2) % 2 == 0 ? sqrt(3) / 2 : 0)
            coords[getidx(x, y)] = (xcoord, ycoord)
        end
    end

    return graph, coords
end
