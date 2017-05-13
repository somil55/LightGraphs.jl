# Betweenness centrality measures
# TODO - weighted, separate unweighted, edge betweenness


@doc_str """
    betweenness_centrality(g[, vs])
    betweenness_centrality(g, k)

Calculate the [betweenness centrality](https://en.wikipedia.org/wiki/Centrality#Betweenness_centrality)
of a graph `g` across all vertices, a specified subset of vertices `vs`, or a random subset of `k`
vertices. Return a vector representing the centrality calculated for each node in `g`.

### Optional Arguments
- `normalize=true`: If true, normalize the betweenness values by the
total number of possible distinct paths between all pairsin the graphs.
For an undirected graph, this number is ``\\frac{(|V|-1)(|V|-2)}{2}``
and for a directed graph, ``{(|V|-1)(|V|-2)}``.
- `endpoints=false`: If true, include endpoints in the shortest path count.

Betweenness centrality is defined as:
``
bc(v) = \\frac{1}{\\mathcal{N}} \sum_{s \\neq t \\neq v}
\\frac{\\sigma_{st}(v)}{\\sigma_{st}}
``.

### References
- Brandes 2001 & Brandes 2008
"""
function parallel_betweenness_centrality(
    g::AbstractGraph,
    vs::AbstractVector = vertices(g);
    normalize=true,
    endpoints=false,verbose=false)

    n_v = nv(g)
    k = length(vs)
    isdir = is_directed(g)


    betweenness = @parallel (+) for s in vs
        if degree(g,s) > 0  # this might be 1?
            state = parallel_dijkstra_shortest_paths(g, s; allpaths=true)
            if endpoints
                parallel_accumulate_endpoints!(state, g, s)
            else
                parallel_accumulate_basic!(state, g, s)
            end
        else
            zeros(Float64,n_v)
        end
    end

    parallel_rescale!(betweenness,
    n_v,
    normalize,
    isdir,
    k)

    return betweenness
end

parallel_betweenness_centrality(g::AbstractGraph, k::Integer; normalize=true, endpoints=false) =
parallel_betweenness_centrality(g, sample(vertices(g), k); normalize=normalize, endpoints=endpoints)



function parallel_accumulate_basic!(
    state::ParallelDijkstraState,
    g::AbstractGraph,
    si::Integer
    )

    n_v = length(state.parents) # this is the ttl number of vertices
    betweenness = zeros(Float64,n_v)
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors

    # make sure the source index has no parents.
    P[si] = []
    # we need to order the source vertices by decreasing distance for this to work.#This is chnged as cmpared to _accumulate_basic!# S = sortperm(state.dists, rev=true)
    S = reverse(state.closest_vertices)
    # S = sortperm(state.dists, rev=true)
    for w in S
        coeff = (1.0 + δ[w]) / σ[w]
        for v in P[w]
            if v > 0
                δ[v] += (σ[v] * coeff)
            end
        end
        if w != si
            betweenness[w] += δ[w]
        end
    end

    return betweenness
end

function parallel_accumulate_endpoints!(
    state::ParallelDijkstraState,
    g::AbstractGraph,
    si::Integer
    )

    n_v = nv(g) # this is the ttl number of vertices
    betweenness = zeros(Float64,n_v)
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors
    v1 = [1:n_v;]
    v2 = state.dists # we need to order the source vertices by decreasing distance for this to work.#This is chnged as cmpared to _accumulate_endpoints!  #
    S = reverse(state.closest_vertices)
    # S = sortperm(state.dists, rev=true)
    s = vertices(g)[si]
    betweenness[s] += length(S) - 1    # 289

    for w in S
        coeff = (1.0 + δ[w]) / σ[w]
        for v in P[w]
            δ[v] += σ[v] * coeff
        end
        if w != si
            betweenness[w] += (δ[w] + 1)
        end
    end

    return betweenness
end

function parallel_rescale!(betweenness::Vector{Float64}, n::Integer, normalize::Bool, directed::Bool, k::Int)
    if normalize
        if n <= 2
            do_scale = false
        else
            do_scale = true
            scale = 1.0 / ((n - 1) * (n - 2))
        end
    else
        if !directed
            do_scale = true
            scale = 1.0 / 2.0
        else
            do_scale = false
        end
    end
    if do_scale
        if k > 0
            scale = scale * n / k
        end
        for v = 1:length(betweenness)
            betweenness[v] *= scale
        end
    end
end
