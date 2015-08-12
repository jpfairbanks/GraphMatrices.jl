module GraphMatrices
if VERSION < v"0.4.0-dev"
	using Docile
	using Compat
end
@doc "A package for using the type system to check types of graph matrices." -> GraphMatrices
import Base.convert
import Base.sparse
import Base.size
import Base.scale
import Base.diag
using FactCheck
import Base.eltype
import Base.size
import Base.ndims
import Base.issym
export  convert,
		SparseMatrix,
		GraphMatrix,
		adjacency,
		Laplacian,
		CombinatorialAdjacency,
		CombinatorialLaplacian,
		NormalizedAdjacency,
		NormalizedLaplacian,
		StochasticAdjacency,
		StochasticLaplacian,
		AveragingAdjacency,
		AveragingLaplacian,
		Noop,
		diag,
		degrees,
		symmetrize,
		prescalefactor,
		postscalefactor



typealias SparseMatrix{T} SparseMatrixCSC{T,Int64}

@doc "An abstract type to allow opertions on any type of graph matrix" ->
abstract GraphMatrix{T}


@doc "The core Adjacency matrix structure. Keeps the vertex degrees around.
Subtypes are used to represent the different normalizations of the adjacency matrix.
Laplacian and its subtypes are used for the different Laplacian matrices.

Adjacency(lapl::Laplacian) provide a generic function for getting the
adjacency matrix of a Laplacian matrix. If your subtype of Laplacian does not provide
an field A for the Adjacency instance, then attach another method to this function to provide
an Adjacency{T} representation of the Laplacian. The Adjacency matrix here 
is the final subtype that corresponds to this type of Laplacian" ->
abstract Adjacency{T} <: GraphMatrix{T}
abstract Laplacian{T} <: GraphMatrix

@doc "Combinatorial Adjacency matrix is the standard adjacency matrix from math" ->
type CombinatorialAdjacency{T,S,V} <: Adjacency{T}
	A::S
	D::V
end

function CombinatorialAdjacency{T}(A::SparseMatrix{T})
	D = vec(sum(A,1))
	return CombinatorialAdjacency{T,SparseMatrix{T},typeof(D)}(A,D)
end


@doc "Normalized Adjacency matrix is \$\\hat{A} = D^{-1/2} A D^{-1/2}\$.
If A is symmetric, then the normalized adjacency is also symmetric
with real eigenvalues bounded by [-1, 1]." ->
type NormalizedAdjacency{T} <: Adjacency{T}
	A::CombinatorialAdjacency{T}
	scalefactor::Vector{T}

	function NormalizedAdjacency(adjmat::CombinatorialAdjacency)
		sf = adjmat.D.^(-1/2)
		return new(adjmat, sf)
	end
end
function NormalizedAdjacency{T}(adjmat::CombinatorialAdjacency{T})
	return NormalizedAdjacency{T}(adjmat)
end

@doc "Transition matrix for the random walk." ->
type StochasticAdjacency{T} <: Adjacency{T}
	A::CombinatorialAdjacency{T}
	scalefactor::Vector{T}

	function StochasticAdjacency(adjmat::CombinatorialAdjacency)
		sf = adjmat.D.^(-1)
		return new(adjmat, sf)
	end
end
function StochasticAdjacency{T}(adjmat::CombinatorialAdjacency{T})
	return StochasticAdjacency{T}(adjmat)
end
@doc "The matrix whos action is to average over each neighborhood." ->
type AveragingAdjacency{T} <: Adjacency{T}
	A::CombinatorialAdjacency{T}
	scalefactor::Vector{T}

	function AveragingAdjacency(adjmat::CombinatorialAdjacency)
		sf = adjmat.D.^(-1)
		return new(adjmat, sf)
	end
end
function AveragingAdjacency{T}(adjmat::CombinatorialAdjacency{T})
	return AveragingAdjacency{T}(adjmat)
end
@doc "Noop: a type to represent don't do anything.
The purpose is to help write more general code for the different scaled GraphMatrix types." ->
type Noop
end

function scale(::Noop, x::Noop)
	return x
end

function scale(::Noop, x::Any)
	return x
end

function scale(x::Any, ::Noop)
	return x
end

function .*(::Noop, x::Any)
	return x
end

function ==(g::GraphMatrix, h::GraphMatrix)
	if typeof(g) != typeof(h)
		return false
	end
	if g.A == h.A
		return true
	end
end

@doc "postscalefactor(M)*M.A*prescalefactor(M) == M " ->
function postscalefactor(::Adjacency)
	return Noop()
end
function postscalefactor(adjmat::NormalizedAdjacency)
	return adjmat.scalefactor
end
function postscalefactor(adjmat::AveragingAdjacency)
	return adjmat.scalefactor
end


@doc "postscalefactor(M)*M.A*prescalefactor(M) == M " ->
function prescalefactor(::Adjacency)
	return Noop()
end
function prescalefactor(adjmat::NormalizedAdjacency)
	return adjmat.scalefactor
end
function prescalefactor(adjmat::StochasticAdjacency)
	return adjmat.scalefactor
end


@doc "Combinatorial Laplacian L = D-A" ->
type CombinatorialLaplacian{T} <: Laplacian{T}
	A::CombinatorialAdjacency{T}
end

@doc "Normalized Laplacian is \$\\hat{L} = I - D^{-1/2} A D^{-1/2}\$.
If A is symmetric, then the normalized Laplacian is also symmetric
with positive eigenvalues bounded by 2." ->
type NormalizedLaplacian{T} <: Laplacian{T}
	A::NormalizedAdjacency{T}
end

@doc "Laplacian version of the StochasticAdjacency matrix." ->
type StochasticLaplacian{T} <: Laplacian{T}
	A::StochasticAdjacency{T}
end

@doc "Laplacian version of the AveragingAdjacency matrix." ->
type AveragingLaplacian{T} <: Laplacian{T}
	A::AveragingAdjacency{T}
end

# function passthrough{T:<GraphMatrix}(f::Function, Type{T})
# 	f(x::Type{T}) = f(x.A)
# end
# eltype(A)	the type of the elements contained in A
# length(A)	the number of elements in A
# ndims(A)	the number of dimensions of A
# size(A)	a tuple containing the dimensions of A
# size(A,n)	the size of A in a particular dimension
# stride(A,k)	the stride (linear index distance between adjacent elements) along dimension k
# strides(A)	a tuple of the strides in each dimension
arrayfunctions = (:eltype, :length, :ndims, :size, :strides, :issym)
for f in arrayfunctions
	@eval $f(a::GraphMatrix) = $f(a.A)
end
size(a::GraphMatrix, i::Integer) = size(a.A, i)
issym(::StochasticAdjacency) = false
issym(::AveragingAdjacency) = false

@doc "degrees of a graph as a Vector." ->
function degrees(adjmat::CombinatorialAdjacency)
	return adjmat.D
end
function degrees(mat::GraphMatrix)
	return degrees(adjacency(mat))
end

# function degrees(lapl::Laplacian)
# 	return degrees(adjacency(lapl))
# end

function adjacency(lapl::Laplacian)
	return lapl.A
end

function adjacency(lapl::GraphMatrix)
	return lapl.A
end

function convert(::Type{Adjacency}, lapl::Laplacian)
	return lapl.A
end
if VERSION < v"0.4"
	combinatorialadjacency(adjmat::Adjacency) = adjmat.A
	# combinatorialadjacency(adjmat::Laplacian) = adjmat.A
end
function convert(::Type{CombinatorialAdjacency}, adjmat::Adjacency)
	return adjmat.A
end

function convert(::Type{SparseMatrix}, adjmat::CombinatorialAdjacency)
	return adjmat.A
end

sparse(m::SparseMatrixCSC{Float64,Int64}) = m

function sparse{M <: Laplacian}(lapl::M)
	adjmat = adjacency(lapl)
	A = sparse(adjmat)
	L = spdiagm(diag(lapl)) - A
	return L
end

function sparse(adjmat::Adjacency)
    A = sparse(adjmat.A)
	return scale(prescalefactor(adjmat), scale(A, postscalefactor(adjmat)))
end

function convert{T}(::Type{SparseMatrix{T}}, adjmat::Adjacency{T})
    A = sparse(adjmat.A)
	return scale(prescalefactor(adjmat), scale(A, postscalefactor(adjmat)))
end

function convert{T}(::Type{SparseMatrix{T}}, lapl::Laplacian{T})
	adjmat = adjacency(lapl)
	A = convert(SparseMatrix{T}, adjmat)
	L = spdiagm(diag(lapl)) - A
	return L
end

function diag(lapl::CombinatorialLaplacian)
	d = lapl.A.D
	return d
end

function diag(lapl::Laplacian)
	return ones(size(lapl)[2])
end


function *{T<:Number}(adjmat::Adjacency{T}, x::Vector{T})
	return  postscalefactor(adjmat) .* (adjmat.A * (prescalefactor(adjmat) .* x))
end

function *{T<:Number}(lapl::Laplacian{T}, x::Vector{T})
	y = adjacency(lapl)*x
	z = diag(lapl) .* x
	return z - y
end

@doc "Symmetrize the matrix.
:triu, :tril, :sum, :or.
use :sum for weighted graphs."->
function symmetrize(A::SparseMatrix, which=:or)
	if which==:triu
		T = triu(A)
	end
	if which==:tril
		T = tril(A)
	end
	if which==:sum
		T = A
	end
	if which==:or
	    M = A + A'
	    M.nzval[M.nzval .== 2] = 1
	else
		M = T + T'
	end
	return M
end

@doc "Only works on Adjacency because the normalizations don't commute with symmetrization"->
function symmetrize(adjmat::CombinatorialAdjacency, which=:or)
	# T = typeof(adjmat)
	# @show T
	# if T <: StochasticAdjacency || T <: AveragingAdjacency
	# 	TypeError("StochasticAdjacency and AveragingAdjacency matrices are nonsymmetric.
	# 		Use NormalizedAdjacency for this purpose")
	# end
	Aprime = symmetrize(adjmat.A, which)
	return CombinatorialAdjacency(Aprime)
end




end
