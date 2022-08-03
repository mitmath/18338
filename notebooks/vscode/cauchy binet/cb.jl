using LinearAlgebra,Combinatorics

n,p = 10,3
A = randn(n,p)
Z = randn(n)
D = Diagonal(1 .+ Z)

𝒮 = combinations(1:n,p) 
𝒫 = powerset(1:n)

sum(  prod(Z[ℐ])*det(A[S,:])^2 for  ℐ∈𝒫, S∈𝒮 if ℐ ⊆ S ) ≈  det(A'*D*A)
