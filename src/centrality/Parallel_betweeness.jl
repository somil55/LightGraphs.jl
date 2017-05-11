
include("/home/divyansh/Desktop/small_graph.jl")

function parallel_betweenness_centrality(
    g::AbstractGraph,
    vs::AbstractVector = vertices(g);
    normalize=true,
    endpoints=false,
    verbose=false)

    vs = vertices(z)
    n_v = nv(z)
    k = length(vs)
    isdir = is_directed(z)
    betweenness = SharedArray{Float64}(n_v)

    @sync @parallel for s in vs
        if degree(z,s) > 0  # this might be 1?
            state = parallel_dijkstra_shortest_paths(z, s; allpaths=true)
            if endpoints
                parallel_accumulate_endpoints!(betweenness, state, z, s)
            else
                parallel_accumulate_basic!(betweenness, state, z, s)
            end
        end
        if verbose
          println(s);
        end
    end

    parallel_rescale!(betweenness,
    n_v,
    normalize,
    isdir,
    k)

    return betweenness
end

@everywhere function parallel_accumulate_basic!(
    betweenness::SharedArray{Float64},
    state::ParallelDijkstraState,
    g::AbstractGraph,
    si::Integer
    )

    n_v = length(state.parents) # this is the ttl number of vertices
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors

    # make sure the source index has no parents.
    P[si] = []

    # we need to order the source vertices by decreasing distance for this to work.
    S = reverse(state.closest_vertices)

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
end


@everywhere function parallel_accumulate_endpoints!(
    betweenness::SharedArray{Float64},
    state::ParallelDijkstraState,
    g::AbstractGraph,
    si::Integer
    )

    n_v = nv(g) # this is the ttl number of vertices
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors
    v1 = [1:n_v;]
    v2 = state.dists
    S = reverse(state.closest_vertices)
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
end


function parallel_rescale!(betweenness::SharedArray{Float64}, n::Integer, normalize::Bool, directed::Bool, k::Int)
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
