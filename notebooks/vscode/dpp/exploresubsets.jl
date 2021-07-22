# explore subsets
using Combinatorics, Printf, LinearAlgebra
N = 5
k = 2

M = randn(N,N)
MI = inv(M)
W = M[:,1:k]
V = MI[1:k,:]'

K = V*W' # (Biorthogonalization says make sure V'W=I, if you haven't already)
V'W


B1 = Matrix( qr(randn(N,k)).Q)
B2 = Matrix( qr(randn(N,k)).Q)
K = B1 * B1' /2 + B2*B2'/2


# Y = randn(N,k)
# Y = Matrix(qr(Y).Q)

# L = Y*Y'
# K = L/(L+I)

 Z = sum( det(K[𝓘,𝓘]) for 𝓘 ∈ combinations(1:N,k))
 Pr(𝓘) = det(K[𝓘,𝓘])

 for j=1:k
    println("---- j=",j)
for A ∈ combinations(1:N,k-1)
#    println(A)    
println(det(K[A,A]), " " ,sum( Pr(𝓘) for 𝓘 ∈ combinations(1:N,k) if A ⊆ 𝓘      ))

#sum( Pr(𝓘) for 𝓘 ∈ combinations(1:N,k) if  [ ]==A ∩ 𝓘   )
 end
end


# for 𝓘∈combinations(1:N,k)
  
#     KI = K[𝓘,𝓘]
#     LI = KI/(I-KI)
#     println("---")
#     println(LI)
#     J = setdiff(1:N,𝓘)
    
#     println(V[J,:]'*W[J,:])
  
# end