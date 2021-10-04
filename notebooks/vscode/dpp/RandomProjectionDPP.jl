using Distributions, Combinatorics

N = 6
n = 3
Y = Matrix(qr(randn(N,n)).Q)
K = Y*Y'

HH(x) = (w = x + norm(x)*[1;zeros(length(x)-1)]; I-2w*w'/w'w)
## the above need not be formed

function randDPP(Y)
    N,n = size(Y)
    𝓘 = Int64[]
    for k=1:n
        p = sum(Y.^2, dims=2)/(n-k+1)
        i = rand( Categorical(p[:]))
        append!(𝓘,[i])   
        Y *= HH( Y[i,:])[:,2:end] 
    end
    return(sort(𝓘))
end

Y = Matrix(qr(randn(N,n)).Q)
K = Y*Y'

## Do a test

hist = Dict( 𝓘=>0 for 𝓘∈combinations(1:N,3) )
t = 100_000
for i=1:t
    hist[randDPP(Y)] += 1
end
hist

for 𝓘∈combinations(1:N,3)
  println(𝓘, " ", hist[𝓘]/t," ",round.(det(K[𝓘,𝓘]),digits=5))
end


