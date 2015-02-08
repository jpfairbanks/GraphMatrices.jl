module TestGraphMatrices
using FactCheck
using GraphMatrices

function subtypepredicate(T)
	pred(x) = issubtype(typeof(x), T)
	return pred
end

function isnot(f::Function)
	return g(x) = !f(x)
end

facts("constructors") do
	n = 10
	mat = SparseMatrix{Float64}(spones(sprand(n,n,0.3)))
	adjmat = CombinatorialAdjacency(mat)
    stochmat = StochasticAdjacency(adjmat)
    adjhat = NormalizedAdjacency(adjmat)
    avgmat = AveragingAdjacency(adjmat)
	context("Adjacency") do
		@fact adjmat.D => vec(sum(mat, 1))
		@fact adjmat.A => mat
		@fact convert(SparseMatrix{Float64}, adjmat) => mat
		@fact convert(SparseMatrix{Float64}, stochmat) => truthy
		@fact convert(SparseMatrix{Float64}, adjhat) => truthy
		@fact convert(SparseMatrix{Float64}, avgmat) => truthy
		@fact stochmat => truthy
		@fact adjhat => truthy
		@fact avgmat => truthy
	end

	context("Laplacian") do
		lapl = CombinatorialLaplacian(CombinatorialAdjacency(mat))
		@fact lapl => truthy
		#constructors that work.
		@fact Adjacency(lapl).A => mat
		@fact StochasticAdjacency(Adjacency(lapl))=> truthy
		@fact NormalizedAdjacency(Adjacency(lapl))=> truthy
		@fact AveragingAdjacency(Adjacency(lapl))=> truthy
		@fact convert(Adjacency, lapl)=> truthy
		@fact convert(SparseMatrix{Float64}, lapl) => truthy
		
		@fact Adjacency(lapl) => subtypepredicate(CombinatorialAdjacency)
		stochlapl = StochasticLaplacian(StochasticAdjacency{Float64}(adjmat))
		@fact Adjacency(stochlapl) => subtypepredicate(StochasticAdjacency)
		averaginglapl = AveragingLaplacian(AveragingAdjacency{Float64}(adjmat))
		@fact Adjacency(averaginglapl) => subtypepredicate(AveragingAdjacency)
		
		normalizedlapl = NormalizedLaplacian(NormalizedAdjacency{Float64}(adjmat))
		@fact Adjacency(normalizedlapl) => subtypepredicate(NormalizedAdjacency)
		@fact Adjacency(normalizedlapl) => isnot(subtypepredicate(CombinatorialAdjacency))

		#constructors that fail.
		@fact_throws CombinatorialAdjacency(lapl)=> truthy
		@fact_throws StochasticLaplacian(lapl) => truthy
		@fact_throws NormalizedLaplacian(lapl) => truthy
		@fact_throws AveragingLaplacian(lapl) => truthy
		@fact_throws convert(CombinatorialAdjacency, lapl)=> truthy
		L = convert(SparseMatrix{Float64}, lapl)
		@fact sum(abs(sum(L,1))) => 0
	end
end

facts("arithmetic") do
	n = 10
	mat = symmetrize(SparseMatrix{Float64}(spones(sprand(n,n,0.3))))
	adjmat = CombinatorialAdjacency(mat)
	stochadj = StochasticAdjacency{Float64}(adjmat)
	averaadj = AveragingAdjacency{Float64}(adjmat)
	hatadj = NormalizedAdjacency(adjmat)
	lapl = CombinatorialLaplacian(adjmat)
	onevec = ones(Float64, n)
	adjmat*ones(Float64, n)
	@fact sum(abs(adjmat*onevec)) => not(0)
	@fact sum(abs(stochadj*onevec/sum(onevec))) => roughly(1)
	@fact sum(abs(lapl*onevec)) => 0
	@pending hatadj => 1
	@pending hatlapl => 1
	@pending averaadj => 1
end

facts("other tests") do
	n = 10
	mat = symmetrize(SparseMatrix{Float64}(spones(sprand(n,n,0.3))))
	adjmat = CombinatorialAdjacency(mat)
	lapl = CombinatorialLaplacian(CombinatorialAdjacency(mat))
	@fact_throws symmetrize(StochasticAdjacency{Float64}(adjmat))
	@fact_throws symmetrize(AveragingAdjacency{Float64}(adjmat))
	@fact_throws symmetrize(NormalizedAdjacency(adjmat)).A => adjmat.A
	@fact symmetrize(adjmat).A => adjmat.A
	
    context("equality testing ") do
        @fact CombinatorialAdjacency(mat) => CombinatorialAdjacency(mat)
        S = StochasticAdjacency(CombinatorialAdjacency(mat))
        @fact S.A => S.A
        @fact Adjacency(S) => not(S.A)
        @fact CombinatorialAdjacency(S) => S.A
        @fact NormalizedAdjacency(adjmat) => not(adjmat)
        @fact StochasticLaplacian(S) => not(adjmat)
        @fact_throws StochasticLaplacian(adjmat) => not(adjmat)
    end
end

end
