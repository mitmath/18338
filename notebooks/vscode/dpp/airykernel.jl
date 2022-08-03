using SpecialFunctions
using Plots
using LinearAlgebra
airy_kernel(x, y) = x==y ? (airyaiprime(x))^2 - x * (airyai(x))^2 :
            (airyai(x) * airyaiprime(y) - airyai(y) * airyaiprime(x)) / (x - y)

s = 3.0
h = .01
𝓘 = 10 : -h : s
𝓘 = rand(7)
## discretize
K = Symmetric(airy_kernel.(𝓘',𝓘))
A = airyai.(𝓘)
A′= airyaiprime.(𝓘)

# Let's try another chocie

A = rand(length(𝓘))
A′ = rand(length(𝓘))
K = (A .* A′' .- A′.*A') ./ (𝓘.-𝓘')
for i=1:size(K,1)
    K[i,i] =0.0
end


ρ = inv(I-K)
R = K/(I-K) #  called the "L" matrix (still divide by det(I+R) to get a prob)
ρ ≈ I+R 

u =  A⋅((I-K)\A) * h
v =  A⋅((I-K)\A′) * h
Q = (I-K)\A 
P = (I-K)\A′
q = Q[end]
p = P[end]  

println( p^2 - s*q^2 - 2*p*q*u + 2*q^2*v)
println(R[end,end])

RR = (Q*P'-P*Q') ./ ( 𝓘.- 𝓘')

