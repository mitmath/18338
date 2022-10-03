using LinearAlgebra, Combinatorics, Distributions, StatsBase

function randprojDPP(Y)

    n = size(Y,2)
    𝓘 = fill(0,n)
    for k=1:n
        p = mean(abs.(Y).^2, dims=2)
        𝓘[k] = rand(Categorical(p[:]))
        Y=(Y*qr(Y[𝓘[k],:]).Q )[:,2:end] 
        #display(Y[𝓘[k],:]) 

    end
    return(sort(𝓘))
end

function randDPP(Y,Λ)
    mask = rand.(Bernoulli.(Λ./(Λ.+1)))
    return(randprojDPP(Y[:,mask]))
end

### test--------------------------------------
N = 2

Λ = rand(N)
Y = Matrix(qr(randn(N,N)).Q)
L = Y * diagm(Λ) * Y'
K = L/(L+I)
randDPP(Y,Λ)

hist = Dict( 𝓘=>0 for 𝓘 ∈ powerset(1:N) )
t = 10000
for i=1:t
    hist[randDPP(Y,Λ)] += 1
end

println("$t trials N=$N")
println("Expnt Theory")
for 𝓘∈powerset(1:N)
   println(round( hist[𝓘]/t,digits=3)," ",round.(det(L[𝓘,𝓘])/det(L+I),digits=3), " ",𝓘)
end

### random wishart query--------------------------------------
# function randwish(N,trials)
# data = fill(0,trials)
# for i=1:trials
#    Y = randn(N,N)+im*randn(N,N)
#    (Λ,X) = eigen(Y*Y')
#    data[i] = length(randDPP(X,Λ))
# end

# c = countmap(data)
# for i ∈ sort( [k for k∈keys(c)])
#     println(i," => ",c[i])
# end

# end

# randwish(100,1000)


# q = rand(5)
# rand.(Bernoulli.(q),10)
