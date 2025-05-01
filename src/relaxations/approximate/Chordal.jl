"""
    chordal_completion!(G::Graphs.SimpleGraph)

Augment `G` by a chordal completion using a greedy minimal fill-in, and also return a perfect elimination ordering.
This is a more efficient implementation of `ChordalGraph.jl/GreedyOrder` (for the "MF" case).

See also [`ChordalGraph.jl`](https://github.com/wangjie212/ChordalGraph)

"""
function chordal_completion!(G::Graphs.SimpleGraph)
    n = Graphs.nv(G)
    H = copy(G)
    order = zeros(Graphs.eltype(G), n)
    @inbounds for i = 1:n
        fill_val = typemax(Int)
        fill_pos = 1
        for j = 1:n
            neib = Graphs.neighbors(H, j)
            lneib = length(neib)
            if lneib > 0
                sg = 0
                for k = 1:lneib, l = k+1:lneib
                    if Graphs.has_edge(H, neib[k], neib[l])
                        sg += 1
                    end
                end
                new_fill_val = ((lneib * (lneib - 1)) >> 1) - sg
                if new_fill_val < fill_val
                    fill_val = new_fill_val
                    fill_pos = j
                end
            elseif order[j] == 0 && fill_val > 0
                fill_val = 0
                fill_pos = j
                break
            end
        end
        order[fill_pos] = i
        neib = copy(Graphs.neighbors(H, fill_pos))
        for j = 1:length(neib)
            Graphs.rem_edge!(H, fill_pos, neib[j])
            for k = j+1:length(neib)
                if Graphs.add_edge!(H, neib[j], neib[k])
                    Graphs.add_edge!(G, neib[j], neib[k])
                end
            end
        end
    end
    return G, order
end

"""
    chordal_cliques!(G::Graphs.SimpleGraph)

Make the given graph chordal, and then calculate its maximal cliques.
This is almost the same implementation as in `ChordalGraph.jl/chordal_cliques!`.

See also [`chordal_completion!`](@ref), [`ChordalGraph.jl`](https://github.com/wangjie212/ChordalGraph)
"""
function chordal_cliques!(G::Graphs.SimpleGraph)
    G, order = chordal_completion!(G)
    n = Graphs.nv(G)
    candidate_cliques = Vector{Vector{Graphs.eltype(G)}}(undef, n)
    for i = 1:n
        candidate_cliques[i] = [intersect(Graphs.neighbors(G, i), findall(x -> order[x] > order[i], 1:n)); i]
        sort!(candidate_cliques[i])
    end
    sort!(candidate_cliques, by=_lensort)
    maximal_cliques = FastVec{eltype(candidate_cliques)}(buffer=length(candidate_cliques))
    for clique in candidate_cliques
        if all(other_clique -> !(clique ⊆ other_clique), maximal_cliques)
            unsafe_push!(maximal_cliques, clique)
        end
    end
    return finish!(maximal_cliques)
end