# GraphMatrices

[![Build Status](https://travis-ci.org/jpfairbanks/GraphMatrices.jl.svg?branch=master)](https://travis-ci.org/jpfairbanks/GraphMatrices.jl)

## Introduction and Motivation

A Julia Package for strongly typed graph matrices.
Many algorithms for analyzing or processing a graph
can be defined in terms of operations on various matrices
associated with the graph. These techniques are typically
found in Spectral Graph Theory.
The Combinatorial Adjacency Matrix has entries A[i,j] = 1 if there is an edge
between vertex i and vertex j.
If D is the matrix where D[i,i] is the degree of vertex i.
Then there are various normalizations to the Adjacency Matrix
that can be used for different purposes.

    StochasticAdjacency = D^{-1}A
    AveragingAdjacency = AD^{-1}
    NormalizedAdjacency = D^{-1/2}AD^{-1/2}

For each of these types there is a corresponding Laplacian matrix.

    CombinatorialLaplacian = D - A
    StochasticLaplacian = I - D^{-1}A
    AveragingLaplacian = I - AD^{-1}
    NormalizedLaplacian = I - D^{-1/2}AD^{-1/2}

The point of this package is to make these available as types.
Because Julia can dispatch on types we can use type information
to make sure that we don't accidentally put a Stochastic Adjacency
into a function that is expecting a Combinatorial Laplacian.

We can also use these designs in order to compute things in a more efficient manner.
The action of each GraphMatrix can be represented in terms of the Adjacency and the normalization.
Thus we can compute the actions without explicitly storing the normalized matrices.
In the future this interface can be extended to include low rank updates to the
adjacency matrix (such as the modularity matrix for community detection).

The eigenvalues  and eigenvectors of the various graph matrices are interdependent
and thus we can compute eigenvalues of one with the most efficient technique
that is available to use. 
For instance the eigenvalues close to 0 for the Normalized Laplacian
have the same eigenvectors as the eigenvalues close to 1 in the Normalized Adjacency Matrix.
Thus by dispatching on types we can compute the small eigenvectors
of the Normalized Laplacian by solving for the large eigenvectors
of the Normalized Adjacency matrix and then transforming the solutions.

## Installation

````julia
    Pkg.add("GraphMatrices")
````

## Types and Usage

This package uses a type hierarchy to represent the 8 differnt graph matrices
all types are parameterized by the type of the entries it contains.

- GraphMatrix
    - Adjacency
        - CombinatorialAdjacency
        - StochasticAdjacency
        - AveragingAdjacency
        - NormalizedAdjacency
    - Laplacian
        - CombinatorialLaplacian
        - StochasticLaplacian
        - AveragingLaplacian
        - NormalizedLaplacian

The combinatorial types can use integer entries but the 3 normalizations require types that support
division or sqrt as appropriate.

### Usage
The usage pattern should be to construct a sparse matrix of some kind usually SparseMatrixCSC{Float64, Int64}
to hold the adjacency matrix data, then wrap it in a CombinatorialAdjacency instance.
Then to explicitly construct the type that is most natural to express your algorithm in.
For example graph partitioning is often expressed in terms of the Fiedler vector which
is the second smallest eigenvector of the Laplacian. So you would implement your algorithm working with the
Laplacian. And define methods that take the Laplacian and convert it to the Adjacency or Normalized Adjacency
when that leads to a faster or more numerically stable implementation.

Math bugs are some of the hardest to detect because the program will run without crashing
but the numbers that it outputs and the conclusion one draws from those numbers will be wrong.
Strong typing should help users of this package catch math bugs as type errors as early as possible.

## The Future!
This package will be simple and provide just a few types and methods in order to make it easier to maintain.
The generalization of this interface to alternative storage formats besides SparseMatrixCSC would be nice
In particular Graphs.jl types could be supported, as could sparse matrix types based on Associative Collections.
This package should be used by any project that is performing any spectal graph algorithms
in order to facilitate development of these important algorithms.
