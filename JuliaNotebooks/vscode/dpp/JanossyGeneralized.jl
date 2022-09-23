N = 7
Y = randn(N,N)
L = Y*Y'
K = L / (L+I)
Pr( 𝓘 ) = det(L[𝓘,𝓘])/det(L+I)

J = [1,2,3,4]
KK = zeros(N,N)
KK[J,J] = K[J,J]
LL=KK/(I-KK)
for ℋ ∈ powerset(J)
  
    println(ℋ," " )
    println(sum( Pr(𝓘) for 𝓘∈powerset(1:N) if  𝓘∩J == ℋ  ))
    println( det(LL[ℋ,ℋ])/det(LL+I))
end