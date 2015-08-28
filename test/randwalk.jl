using GraphMatrices
using FactCheck
if VERSION < v"0.4.0-dev"
	using Docile
	using Compat
end

@doc "Computes the stationary distribution of a random walk" ->
function stationarydistribution(R::StochasticAdjacency; kwargs...)
	er = eigs(R, nev=1, which=:LR; kwargs...)
	l1 = er[1][1]
	abs(l1 -1) < 1e-8 || error("failed to compute stationary distribution")
	p = real(er[2][:,1])
	if p[1] < 0
		for i in 1:length(p)
			p[i] = -p[i]
		end
	end
	return p
end

function stationarydistribution(A::CombinatorialAdjacency; kwargs...)
	R = StochasticAdjacency(A)
	stationarydistribution(R; kwargs...)
end

facts("Random Walk Demo") do
	n = 100
	p = 16/n
	M = sprand(n,n, p)
	M.nzval[:] = 1.0
	A = CombinatorialAdjacency(M)
	@fact stationarydistribution(A; ncv=10) --> x->(all(x.>=0))
end
