<<<<<<< HEAD
using LinearAlgebra, Statistics, Combinatorics, Distributions

HH(x) = (w = x + norm(x)*[1;zeros(length(x)-1)]; I-2w*w'/w'w)
=======
using LinearAlgebra, Combinatorics, Distributions, StatsBase
>>>>>>> df5b3b8d70a2d067c5e14563c9190335b937ecd1

function randprojDPP(Y)
    n = size(Y,2)
    𝓘 = fill(0,n)
    for k=1:n
        p = mean(abs.(Y).^2, dims=2)
        𝓘[k] = rand(Categorical(p[:]))
        Y=(Y*qr(Y[𝓘[k],:]).Q )[:,2:end] # just a householder matrix 
    end
    return(sort(𝓘))
end

function randDPP(Y,Λ)
    mask = rand.(Bernoulli.(Λ./(Λ.+1)))
    return(randprojDPP(Y[:,mask]))
end

### test--------------------------------------
# N = 5
# Λ = rand(N)
# Y = Matrix(qr(randn(N,N)).Q)
# L = Y * diagm(Λ) * Y'

# hist = Dict( 𝓘=>0 for 𝓘 ∈ powerset(1:N) )
# t = 500_000
# for i=1:t
#     hist[randDPP(Y,Λ)] += 1
# end

# for 𝓘∈powerset(1:N)
#    println(round( hist[𝓘]/t,digits=3)," ",round.(det(L[𝓘,𝓘])/det(L+I),digits=3), " ",𝓘)
# end

### random wishart query--------------------------------------
function randwish(N,trials)
data = fill(0,trials)
for i=1:trials
   Y = randn(N,N)+im*randn(N,N)
   (Λ,X) = eigen(Y*Y')
   data[i] = length(randDPP(X,Λ))
end

c = countmap(data)
for i ∈ sort( [k for k∈keys(c)])
    println(i," => ",c[i])
end

end

randwish(100,10)

