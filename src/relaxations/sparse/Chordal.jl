"""
    chordal_cliques!(graphs::Graphs.SimpleGraph; alg::CliqueTrees.EliminationAlgorithm=CliqueTrees.MF())

Make the given graph chordal, and then calculate its maximal cliques.
"""
function chordal_cliques!(graph::Graphs.SimpleGraph; alg::CliqueTrees.EliminationAlgorithm=CliqueTrees.MF())
    # compute a tree decomposition using the given algorithm
    order, tree = CliqueTrees.cliquetree(graph; alg)

    # compute maximal cliques
    cliques = Vector{Vector{Int}}(undef, length(tree))

    for (i, clique) in enumerate(tree)
        # do these need to be sorted?
        cliques[i] = sort!(order[clique])
    end

    # triangulate graph
    filledgraph = CliqueTrees.FilledGraph(tree)

    for v in Graphs.vertices(filledgraph), w in Graphs.neighbors(filledgraph, v)
        Graphs.add_edge!(graph, order[v], order[w])
    end

    return cliques
end

"""
    chordal_cliques(graph::Graphs.SimpleGraph; alg::CliqueTrees.EliminationAlgorithm=CliqueTrees.MF())

Non-mutating version of [`chordal_cliques!`](@ref)
"""
function chordal_cliques(graph::Graphs.SimpleGraph; alg::CliqueTrees.EliminationAlgorithm=CliqueTrees.MF())
    # compute a tree decomposition using the given algorithm
    order, tree = CliqueTrees.cliquetree(graph; alg)

    # compute maximal cliques
    cliques = Vector{Vector{Int}}(undef, length(tree))

    for (i, clique) in enumerate(tree)
        # do these need to be sorted?
        cliques[i] = sort!(order[clique])
    end

    return cliques
end
