#Stores the parallel Heap entry
struct ParallelDijkstraHeapEntry{T, U<:Integer}
    vertex::U
    dist::T
end

isless(e1::ParallelDijkstraHeapEntry, e2::ParallelDijkstraHeapEntry) = e1.dist < e2.dist

"""
    struct DijkstraState{T, U}

An [`AbstractPathState`](@ref) designed for Dijkstra shortest-paths calculations.
"""

abstract type ParallelAbstractPathState end

struct ParallelDijkstraState{T, U<:Integer}<: ParallelAbstractPathState
    parents::Vector{U}
    dists::Vector{T}
    predecessors::Vector{Vector{U}}
    pathcounts::Vector{U}
    closest_vertices::Vector{U}
end

"""
    dijkstra_shortest_paths(g, srcs, distmx=DefaultDistance());

Perform [Dijkstra's algorithm](http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm)
on a graph, computing shortest distances between `srcs` and all other vertices.
Return a [`DijkstraState`](@ref) that contains various traversal information.

### Optional Arguments
- `allpaths=false`: If true, returns a [`DijkstraState`](@ref) that keeps track of all
predecessors of a given vertex.
"""
function parallel_dijkstra_shortest_paths(
    g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T}=DefaultDistance();
    allpaths=false
    ) where T where U<:Integer
    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    closest_vertices = Vector{U}()
    preds = fill(Vector{U}(),nvg)
    visited = zeros(Bool, nvg)
    pathcounts = zeros(Int, nvg)
    H = Vector{ParallelDijkstraHeapEntry{T, U}}()  # this should be Vector{T}() in 0.4, I think.
    dists[srcs] = zero(T)
    pathcounts[srcs] = 1

    sizehint!(H, nvg)
    sizehint!(closest_vertices, nvg)

    for v in srcs
        heappush!(H, ParallelDijkstraHeapEntry{T, U}(v, dists[v]))
        visited[v] = true
    end

    while !isempty(H)
        hentry = heappop!(H)
            # info("Popped H - got $(hentry.vertex)")
        u = hentry.vertex
        push!(closest_vertices, u)
        for v in out_neighbors(g,u)
            alt = (dists[u] == typemax(T))? typemax(T) : dists[u] + distmx[u,v]
            if !visited[v]
                dists[v] = alt
                parents[v] = u
                pathcounts[v] += pathcounts[u]
                visited[v] = true
                if allpaths
                    preds[v] = [u;]
                end
                heappush!(H, ParallelDijkstraHeapEntry{T, U}(v, alt))
                # info("Pushed $v")
            else
                if alt < dists[v]
                    dists[v] = alt
                    pathcounts[v] = 0
                    preds[v] = []
                    parents[v] = u
                    heappush!(H, ParallelDijkstraHeapEntry{T, U}(v, alt))
                end
                if alt == dists[v]
                    pathcounts[v] += pathcounts[u]
                    if allpaths
                        push!(preds[v], u)
                    end
                end
            end
        end
    end


    for s in vertices(g)
      if !visited[s]
        push!(closest_vertices,s)
      end
    end

    pathcounts[srcs] = 1
    parents[srcs] = 0
    for src in srcs
        preds[src] = []
    end

    return ParallelDijkstraState{T, U}(parents, dists, preds, pathcounts,closest_vertices)
end

parallel_dijkstra_shortest_paths(g::AbstractGraph, src::Integer, distmx::AbstractMatrix = DefaultDistance(); allpaths=false) =
parallel_dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths)
