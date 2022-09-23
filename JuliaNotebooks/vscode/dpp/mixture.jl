using StatsBase
import Base.rand

struct Mixture
    M :: Matrix
    p :: Vector
end

function rand(m::Mixture)
   j = sample( Weights(m.p) )
   i = sample( Weights(m.M[:,j])) 
end

n = 5
M = rand(n,n); M ./= sum(M,dims=1)
p = rand(n); p /= sum(p)
t = 1e6
samples = [rand(Mixture(M,p)) for k=1:t]

theory = M*p
experiment = [ mean(samples.==i) for i=1:n]

[theory experiment]