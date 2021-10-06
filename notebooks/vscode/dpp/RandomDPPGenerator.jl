using LinearAlgebra, Statistics, Combinatorics, Distributions

HH(x) = (w = x + norm(x)*[1;zeros(length(x)-1)]; I-2w*w'/w'w)

function randprojDPP(Y)
    n = size(Y,2)
    𝓘 = fill(0,n)
    for k=1:n
        p = sum(Y.^2, dims=2)/(n-k+1)
        𝓘[k] = rand( Categorical(p[:]))
        Y *=  HH( Y[𝓘[k],:])[:,2:end] 
    end
    return(sort(𝓘))
end

function randDPP(Y,Λ)
    mask = rand.(Bernoulli.(Λ./(Λ.+1)))
    return(randprojDPP(Y[:,mask]))
end

### test--------------------------------------
N = 4
Λ = rand(N)
Y = Matrix(qr(randn(N,N)).Q)
L = Y * diagm(Λ) * Y'

hist = Dict( 𝓘=>0 for 𝓘∈powerset(1:N) )
t = 100_000
for i=1:t
    hist[randDPP(Y,Λ)] += 1
end
hist
for 𝓘∈powerset(1:N)
  println(𝓘, " ", hist[𝓘]/t," ",round.(det(L[𝓘,𝓘])/det(L+I),digits=5))
end
