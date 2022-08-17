
export Coord, get_neighbor_graph, nnodes, build_graph
export GraphParams, ChainGraphParams, TriangularGraphParams, SquareGraphParams, HoneycombGraphParams

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
############################################################################################

abstract type GraphParams end

@with_kw struct ChainGraphParams <: GraphParams
    N::Int
    periodic::Bool=false
end

@with_kw struct TriangularGraphParams <: GraphParams
    Nx::Int
    Ny::Int
    yperiodic::Bool=false
end

@with_kw struct HoneycombGraphParams <: GraphParams
    Nx::Int
    Ny::Int
    yperiodic::Bool=false
end

@with_kw struct SquareGraphParams <: GraphParams
    Nx::Int
    Ny::Int
    yperiodic::Bool=false
end

# Get number of sites
nnodes(params::GraphParams) = @error "nnodes for type $(typeof(params)) is not implemented"
nnodes(params::ChainGraphParams) = params.N
nnodes(params::TriangularGraphParams) = params.Nx * params.Ny
nnodes(params::SquareGraphParams) = params.Nx * params.Ny
nnodes(params::HoneycombGraphParams) = params.Nx * params.Ny

# Build Graphs

build_graph(params::GraphParams) = @error "build_graph for type $(typeof(params)) is not implemented"

function build_graph(params::ChainGraphParams)
    @unpack N, periodic = params

    coords = [(i - 1, 0) for i in 1:N]

    edges = Edge.([(i, i+1) for i in 1:(N-1)])
    periodic && push!(edges, Edge(N, 1))

    graph = SimpleGraph(edges)

    return graph, coords
end

function build_graph(params::TriangularGraphParams)
    @unpack Nx, Ny, yperiodic = params

    yperiodic = yperiodic && (Ny > 2)

    n = Nx * Ny
    nbonds = 3n - 2Ny + (yperiodic ? 0 : -2Nx + 1)
    edges = Vector{Tuple{Int, Int}}(undef, nbonds)
    coords = Vector{Coord}(undef, n)

    b = 0
    for i in 1:n
        x = (i - 1) ÷ Ny + 1
        y = mod1(i, Ny)

        # set coordinates
        coords[i] = (x - 1, y - 1)

        # x-direction bonds
        if x < Nx
            edges[b += 1] = (i, i + Ny)
        end

        # 2d bonds
        if Ny > 1
            # vertical / y-periodic diagonal bond
            if (i + 1 <= n) && ((y < Ny) || yperiodic)
                edges[b += 1] = (i, i + 1)
            end
            # periodic vertical bond
            if yperiodic && y == 1
                edges[b += 1] = (i, i + Ny - 1)
            end
            # diagonal bonds
            if x < Nx && y < Ny
                edges[b += 1] = (i, i + Ny + 1)
            end
        end
    end

    graph = SimpleGraph(Edge.(edges))

    return graph, coords
end

function build_graph(params::SquareGraphParams)
    @unpack Nx, Ny, yperiodic = params

    yperiodic = yperiodic && (Ny > 2)

    n = Nx * Ny
    nbonds = 2n - Ny + (yperiodic ? 0 : -Nx)
    edges = Vector{Tuple{Int, Int}}(undef, nbonds)
    coords = Vector{Coord}(undef, n)

    b = 0
    for i in 1:n
        x = (i - 1) ÷ Ny + 1
        y = mod1(i, Ny)

        # set coordinates
        coords[i] = (x - 1, y - 1)

        if x < Nx
            edges[b += 1] = (i, i + Ny)
        end

        if Ny > 1
            if y < Ny
                edges[b += 1] = (i, i + 1)
            elseif yperiodic
                edges[b += 1] = (i, i - Ny + 1)
            end
        end
    end

    graph = SimpleGraph(Edge.(edges))

    return graph, coords
end

"""
Build honeycomb lattice in armchair configuration
"""
function build_graph(params::HoneycombGraphParams)
    @unpack Nx, Ny, yperiodic = params

    getidx(x, y) = (x - 1) * Ny + y
    wrap(x, k) = mod1(mod1(x, k) + k, k)

    # Build graph
    yperiodic = yperiodic && (Ny > 2)

    edges = Vector{Tuple{Int, Int}}()

    function add_edge!(x1, y1, x2, y2)
        idx1 = getidx(x1, y1)
        idx2 = getidx(x2, y2)
        push!(edges, (idx1, idx2))
    end

    for x in 1:(Nx - 1)
        for y in 1:Ny
            # bind to next row
            add_edge!(x, y, x + 1, y)

            # Alternate between linking up and linking down on odd rows
            if x % 2 == 1
                crossy = y + ((x ÷ 2) % 2 == 0 ? 1 : -1)
                if crossy ∈ 1:Ny
                    add_edge!(x, y, x + 1, crossy)
                elseif yperiodic
                    add_edge!(x, y, x + 1, wrap(crossy, Ny))
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

    coords = Vector{Coord}(undef, Nx * Ny)
    for x in 1:Nx
        for y in 1:Ny
            xcoord = (x - 1) * 3/4 + ((x % 2) == 0 ? -1/4 : 0)
            ycoord = (y - 1) * sqrt(3) + ((x ÷ 2) % 2 == 0 ? sqrt(3) / 2 : 0)
            coords[getidx(x, y)] = (xcoord, ycoord)
        end
    end

    return graph, coords
end
