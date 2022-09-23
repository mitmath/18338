using Distributions, StatsBase

# λ1 = 3
# λ2 = 1

# {1} with prob .75
# {2} with prob .4
1+1

I = []
if rand(Bernoulli(.75))  append!(I,1) end
if rand(Bernoulli(.5))  append!(I,2) end
I

# pr []  = .25 * .6
# pr 1  = .75 * .6
# pr 2 = .25 * .4
# pr 1,2 = .75 * .4