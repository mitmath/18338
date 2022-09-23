# "broadcasting  - the "." operator, the point operator, pointwise or elementwise operations"

a = rand(3,3)
sin.(a)
exp.(a) # exp(a) is the "matrix exponential"

i = 1:4
i' .+ i
i' .* i

typeof.([1,1.0,"one",ones(Int,3,3)])

a = rand(4,4)
b = rand(1,4)
#size( a .+ b )
a = [1 2 3; 4 5 6; 7 8 9]
b = [1 10 100]
