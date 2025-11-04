using LinearAlgebra, Combinatorics, Distributions, StatsBase

function randprojDPP(Y)

    n = size(Y,2)
    𝓘 = fill(0,n)
    for k=1:n
        p = mean(abs.(Y).^2, dims=2)
        𝓘[k] = rand(Categorical(p[:]))
        Y=(Y*qr(Y[𝓘[k],:]).Q )[:,2:end] 
    end
    return(sort(𝓘))
end

function randDPP(Λ,Q)
    mask = rand.(Bernoulli.(Λ./(Λ.+1)))
    return(randprojDPP(Q[:,mask]))
end

randDPP(L::Matrix) = randDPP(eigen(L)...)


