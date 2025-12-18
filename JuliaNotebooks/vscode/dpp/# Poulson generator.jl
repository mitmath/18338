# Poulson generator
using Distributions, LinearAlgebra, Combinatorics

function randDPP_P(Kernel)
   K = copy(Kernel)
   N = size(K,2) 
   𝓘 =  Int[]
   for i=1:N
     rand(Bernoulli(K[i,i])) ?  push!(𝓘 ,i) :  K[i,i] -= 1
     K[i+1:N,i+1:N] -= K[i,i+1:N]*K[i+1:N,i]'/K[i,i]
   end
   return 𝓘 
end

### test--------------------------------------
N = 3
L = randn(N,N)
L *= L'
K = L/(L+I)


hist = Dict( 𝓘=>0 for 𝓘 ∈ powerset(1:N) )
t = 10000
for i=1:t
    hist[randDPP_P(K)] += 1
end

println("$t trials N=$N")
println("Expmnt Theory")

pp(𝓘) = length(𝓘)>0 ? 𝓘 : "∅" # pretty print the empty set

for 𝓘∈powerset(1:N)
   println(round( hist[𝓘]/t,digits=3)," ",round.(det(L[𝓘,𝓘])/det(L+I),digits=3), " ",pp(𝓘))
end

