using Combinatorics, LinearAlgebra, NamedArrays
# worth seeing Rota: https://www.maths.ed.ac.uk/~v1ranick/papers/rota1.pdf
# or Stanley: http://www-math.mit.edu/~rstan/ec/ec1.pdf#page=303
# also wikipedia: https://en.wikipedia.org/wiki/M%C3%B6bius_inversion_formula#On_posets
# and wikipedia on inclusion-exclusion: https://en.wikipedia.org/wiki/Inclusion%E2%80%93exclusion_principle#Other_forms

𝒫(N) = collect(powerset(1:N))
ζ(N) = Int.([ A ⊆ B for A∈𝒫(N), B∈𝒫(N) ])
μ(N) =   Int.(inv(ζ(N)))

N = 3
Y = randn(N,N)
L = Y'Y
K = L/(L+I)

PDF = [det(L[𝓘,𝓘])/det(L+I) for 𝓘∈𝒫(N)]
CDF = [det(K[𝓘,𝓘])  for 𝓘∈𝒫(N)]
(CDF  ≈ ζ(N) * PDF)  && (PDF  ≈ μ(N) * CDF)

# labels = string.(𝒫(N));labels[1] = "∅"
# display(NamedArray( ζ(N) , (labels,labels) ) )
# display(NamedArray( μ(N) , (labels,labels) ) )