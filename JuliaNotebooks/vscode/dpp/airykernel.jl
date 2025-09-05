#using SpecialFunctions
#using Plots
using LinearAlgebra
# airy_kernel(x, y) = x==y ? (airyaiprime(x))^2 - x * (airyai(x))^2 :
#             (airyai(x) * airyaiprime(y) - airyai(y) * airyaiprime(x)) / (x - y)

# s = 3.0
# h = .01
# 𝓘 = 10 : -h : s
 
# ## discretize
# K = Symmetric(airy_kernel.(𝓘',𝓘))
# A = airyai.(𝓘)
# A′= airyaiprime.(𝓘)

#𝓘 = rand(3)
#A = randn(length(𝓘))
#A′ = randn(length(𝓘))



n = 3
𝓘 = rand(n)
A = rand(n)
A′ =  rand(n)
K = (A .* A′' .- A′.*A') ./ (𝓘.-𝓘')
for i=1:length(A)
    K[i,i] =0.0
end
R = simplify.(K/(I-K))
Q = (I-K)\A 
P = (I-K)\A′
#q = Q[end]
#p = P[end]
RR = simplify.((Q*P'-P*Q') ./ ( 𝓘.- 𝓘'))
display(RR)
display(R)
E = Matrix(I,n,n)
Δ = kron(E,Diagonal(𝓘)) .- kron( Diagonal(𝓘),E)


#ρ = inv(I-K)
 #  called the "L" matrix (still divide by det(I+R) to get a prob)
#ρ ≈ I+R 

#u =  A⋅((I-K)\A) * h
#v =  A⋅((I-K)\A′) * h
  

#println( p^2 - s*q^2 - 2*p*q*u + 2*q^2*v)
#println(R[end,end])


