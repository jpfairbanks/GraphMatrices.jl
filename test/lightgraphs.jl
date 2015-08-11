using LightGraphs
using GraphMatrices
using FactCheck

import GraphMatrices.CombinatorialAdjacency

function CombinatorialAdjacency(A)
    D = indegree(A, vertices(A))
    return CombinatorialAdjacency{Float64, typeof(A), typeof(D)}(A,D)
end

facts("constructors") do
    g = PathGraph(10)
    Ag = CombinatorialAdjacency(g)
    Am = CombinatorialAdjacency(adjacency_matrix(g))
    @fact Ag => truthy
end
