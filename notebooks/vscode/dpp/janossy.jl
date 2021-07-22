#Janossy Densities: https://arxiv.org/abs/math-ph/0212063
# Formula (6) on p.3

# The key idea is to start with a projection DPP 
# which means the probability space consists of sets of size n
# then we intersect with a specfied set 𝓘Janossy
# so now we have a probality space on the powerset of 𝓘Janossy
# and the result is that this is a DPP 

# Generate a Random Projection DPP
using Combinatorics, Printf, LinearAlgebra
N = 10
n = 4
Y = Matrix(qr(randn(N,n)).Q)

K = Y*Y'
L = zero.(K)
𝓘Janossy = [1,3,4] # Restrict to an "I" 

L[𝓘Janossy,𝓘Janossy] = K[𝓘Janossy,𝓘Janossy]
L = L/(I-L) # Formula (6) on p.3 

println("====")
println("Janossy = ", 𝓘Janossy)
probs=Float64[]
for 𝓘 ∈ powerset(𝓘Janossy)
    println("𝓘 = ",𝓘)
    p = sum( Float64[det(K[J,J]) for J ∈ combinations(1:N,n)  if  J ∩  𝓘Janossy == 𝓘] ) # determinants of size n
    q = det(L[𝓘,𝓘])/det(I+L)
    println(p," ",q)
    push!(probs,p)
end  
println("total probability = ",sum(probs))
