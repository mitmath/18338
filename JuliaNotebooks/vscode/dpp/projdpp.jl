using Combinatorics, LinearAlgebra
N,n = 6,3
Y,_ = qr(randn(N,n))
Y = Matrix(Y)
L = Y*Y'
# This matrix L is a projection


d = Dict{Array{Int64,1},Float64}()
 for 𝓘 ∈ combinations(1:N, n)
    prob = det(L[𝓘,𝓘]) / binomial(N,n)
    for 𝒬 ∈ combinations(𝓘)
        if haskey(d,𝒬)
            d[𝒬] += prob
        else
            d[𝒬] = prob
        end
    end
end

d = Dict{Array{Int64,1},Float64}()
 for 𝓘 ∈ combinations(1:N, n)
    prob = det(L[𝓘,𝓘]) / binomial(N,n)
    for 𝒬 ∈ combinations(𝓘)
        if haskey(d,𝒬)
            d[𝒬] += prob
        else
            d[𝒬] = prob
        end
    end
end

for k=1:n
for 𝒬 ∈ combinations(1:N,k)    
    println(𝒬," ",det(L[𝒬,𝒬])/binomial(N,n) ," ", d[𝒬]  ) 
end
end


[[sum(det(L[𝒬,𝒬]) for 𝒬 ∈ combinations(1:N,k) ) for k=0:n] [binomial(n,k) for k=0:n]]