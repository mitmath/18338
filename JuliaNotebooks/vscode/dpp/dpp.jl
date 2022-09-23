using LinearAlgebra, Combinatorics, Printf
N = 4

A = randn(N,N)
B = randn(N,N)

function mixAandB(A,B,𝓘)
    C = copy(A)
    C[:,𝓘] .= B[:,𝓘]  
    return C
end

## Power set
𝒫 = powerset(1:N)


## Check identity
det(A+B), sum( det(mixAandB(A,B,𝓘)) for 𝓘 ∈ 𝒫 )


## Generate a random DPP
Y = randn(N,N)
L = Y'Y
for 𝓘 ∈ 𝒫
    @printf("  %10s :", 𝓘) 
    @printf(" %s\n",  det(L[𝓘,𝓘])/det(L+I) )
end
𝒫 = append!([[]],[A for A in combinations(1:N)])
𝒫 =[𝒫;]
print("         Sum : ")
sum(  det(L[𝓘,𝓘])/det(L+I)  for 𝓘 ∈ 𝒫 )

H = [2,4]
Hc = setdiff(1:N,H)
println( sum(  det(L[𝓘,𝓘])/det(L+I)  for 𝓘 ∈ 𝒫  if 𝓘 ⊆ Hc) )
K = L/(I+L)
IH = Diagonal(  1:N .∈ [H])
IHc = I -IH
println( det(I-K[H,H]))
println( det(I+IHc*L)/det(I+L))
# println("-------")
# MU = Integer.([ A ⊆ B for A∈𝒫, B∈𝒫 ])
# display(MU)
# ZETA = Integer.(inv(MU))