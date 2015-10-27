"""Performs a random walk on graph `g` starting at vertex `s` and continuing for
a maximum of `niter` steps. Returns a vector of vertices visited in order.
"""
function randomwalk(g::SimpleGraph, s::Integer, niter::Integer)
s in vertices(g) || throw(BoundsError())
visited = Vector{Int}()
sizehint!(visited, niter)
currs = s
i = 1
while i <= niter
    push!(visited, currs)
    i += 1
    nbrs = out_neighbors(g,currs)
    length(nbrs) == 0 && break
    currs = rand(nbrs)
end
return visited[1:i-1]
end

"""Performs a [self-avoiding walk](https://en.wikipedia.org/wiki/Self-avoiding_walk)
on graph `g` starting at vertex `s` and continuing for a maximum of `niter` steps.
Returns a vector of vertices visited in order.
"""
function saw(g::SimpleGraph, s::Integer, niter::Integer)
s in vertices(g) || throw(BoundsError())
visited = Vector{Int}()
svisited = Set{Int}()
sizehint!(visited, niter)
sizehint!(svisited, niter)
currs = s
i = 1
while i <= niter
    push!(visited, currs)
    push!(svisited, currs)
    i += 1
    nbrs = setdiff(Set(out_neighbors(g,currs)),svisited)
    length(nbrs) == 0 && break
    currs = rand(collect(nbrs))
end
return visited[1:i-1]
end

"""Returns the graph bloom of a graph `g` starting with vertex `s` and continuing
for `niter` iterations. Returns a vector of vectors for vertices visited at each
iteration."""
function bloom(g::SimpleGraph, s::Integer, niter::Integer)
    bloomvec = Vector{Vector{Int}}()
    itervec = Vector{Int}()
    sizehint!(bloomvec, niter)
    # for i = 1:niter
    #     push!(bloomvec, Vector{Int}())
    # end
    push!(bloomvec, [s])
    for i = 2 : niter
        newneighbors = unique(vcat(bloomvec[i-1], [out_neighbors(g,x) for x in bloomvec[i-1]]...))
        length(newneighbors) == length(bloomvec[i-1]) && break
        push!(bloomvec,newneighbors)
    end
    return bloomvec
end

# recursive bloom at level iternum: returns srcs for iternum+1.
function bloom!(
    nvec::Vector{Vector{Int}},
    g::SimpleGraph,
    srcs::Vector{Int},
    itervec::Int,
    itermax::Int
    )
    itervec > itermax && return nothing
    println("srcs = $srcs, nvec[itervec-1] = $(nvec[itervec-1])")
    (Set(srcs) == Set(nvec[itervec-1])) && return nothing
    nvec[itervec] = union(nvec[itervec-1],srcs)
    sset = Set{Int}()
    for i in srcs
        union!(sset, out_neighbors(g,i))
    end
    bloom!(nvec, g, collect(sset), itervec+1, itermax)
    return nothing
end
