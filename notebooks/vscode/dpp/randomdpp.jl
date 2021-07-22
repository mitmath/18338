# Generate a random DPP
using Combinatorics, Printf, LinearAlgebra
N = 5
Y = randn(N,N)
L = Y'Y
K = L/(L+I)
for 𝓘 ∈ powerset(1:N)
    @printf("  %10s :", 𝓘) 
    @printf(" %s\n",  det(L[𝓘,𝓘])/det(L+I) )
end

  print("         Sum : ")
  println(sum(  det(L[𝓘,𝓘])/det(L+I)  for 𝓘 ∈ powerset(1:N) ))


for 𝓘 ∈ powerset(1:N)
    v = zeros(N)
    v[𝓘] .= 1
     I𝓘 = Diagonal(  1 .- v     )
    println(  (-1)^length(𝓘)  * det(I𝓘 - K), " ",det(L[𝓘,𝓘])/det(I+L))
 end