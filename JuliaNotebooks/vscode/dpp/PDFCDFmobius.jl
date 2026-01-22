using Combinatorics, LinearAlgebra, NamedArrays

𝒫(N) = collect(powerset(1:N))

N = 3
Y = randn(N,N)
L = Y'Y
K = L/(L+I)
Id(𝓘) = Diagonal([i∈𝓘 for i=1:N])
Īd(𝓘) = I - Id(𝓘)

PDF1 = [det(L[𝓘,𝓘])/det(L+I) for 𝓘∈𝒫(N)]
PDF2 = [ (-1)^length(𝓘) * det( Īd(𝓘) - K) for 𝓘∈𝒫(N)]
PDF3 = [ sum( (-1)^length(setdiff(ℑ,𝓘)) *det(K[ℑ,ℑ])  for ℑ∈𝒫(N) if ℑ ⊇ 𝓘 ) for 𝓘∈𝒫(N)   ]
[PDF1 PDF2 PDF3]
#ζ(N) = Int.([ A ⊆ B for A∈𝒫(N), B∈𝒫(N) ])
#μ(N) =   Int.(inv(ζ(N)))
# worth seeing Rota: https://www.maths.ed.ac.uk/~v1ranick/papers/rota1.pdf
# or Stanley: http://www-math.mit.edu/~rstan/ec/ec1.pdf#page=303
# also wikipedia: https://en.wikipedia.org/wiki/M%C3%B6bius_inversion_formula#On_posets
# and wikipedia on inclusion-exclusion: https://en.wikipedia.org/wiki/Inclusion%E2%80%93exclusion_principle#Other_forms

# labels = string.(𝒫(N));labels[1] = "∅"
# display(NamedArray( ζ(N) , (labels,labels) ) )
# display(NamedArray( μ(N) , (labels,labels) ) )
# CDF = [det(K[𝓘,𝓘])  for 𝓘∈𝒫(N)]
# (CDF  ≈ ζ(N) * PDF1)  && (PDF1  ≈ μ(N) * CDF)
