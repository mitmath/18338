# Generate a random DPP
using Combinatorics, Printf, LinearAlgebra
N = 2
Y = randn(N,N) # not orthogonal
L = Y'Y # not projection
K = L/(L+I)
  
# pdfs for 𝓘
for 𝓘 ∈ powerset(1:N)
    @printf("  %10s :", 𝓘) 
    @printf(" %s ",  det(L[𝓘,𝓘])/det(L+I) )
    𝓘c = setdiff( 1:N,𝓘 )
    @printf(" %s\n",  det(inv(L)[𝓘c,𝓘c])/det(inv(L)+I) )  # complementary DPP
end

  print("         Sum : ")
  println(sum(  det(L[𝓘,𝓘])/det(L+I)  for 𝓘 ∈ powerset(1:N) ))

Pr(𝓘) = det(L[𝓘,𝓘])/det(L+I) 

# cdfs for 𝓘
println("cdfs")
for H ∈ powerset(1:N)
   print( sum( Pr(𝓘) for 𝓘∈ powerset(1:N) if 𝓘 ⊇  H) , " ")
   print( det(K[H,H]), " ")
   Hc = setdiff(1:N, H)
   
   println( sum( Pr(setdiff(1:N,𝓘)) for 𝓘∈ powerset(1:N) if 𝓘 ⊆ setdiff(1:N, H)) , " ")
 end


# for H ∈ powerset(1:N)
#   print( sum( Pr(𝓘) for 𝓘∈ powerset(1:N) if 𝓘∩ H  == []) , " ")
#   II = setdiff(1:N,H)
#   println( det(I-K[H,H])) 
# end
println()

for 𝓘 ∈ powerset(1:N)

  print( sum( det(L[J,J])/det(L+I) for J∈ powerset(1:N) if J ⊇ 𝓘), " ")
  print( det(K[𝓘,𝓘]), " ")
  print(  sum( det(L[J,J])/det(L+I) for J∈ powerset(1:N) if J ⊆ 𝓘), " ") 
  𝓘c = setdiff(1:N,𝓘)
  print( det( (I-K)[𝓘c,𝓘c]))
  println()
end