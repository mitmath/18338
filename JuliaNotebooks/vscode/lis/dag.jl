using LinearAlgebra, Plots, SparseArrays, AbstractAlgebra

p = [3,1,4,5,9,2,6,8,7]
p = [8,3,7,9,2,5,4,1,10,6]
A = triu(p'.>p)

λ = Partition(collect(12:-1:1))

y = YoungTableau([5,4,2,2,1])

matrix_repr(y)
m = Matrix(y)



y = YoungTableau([3,2,1])
YoungTableau(Partition([3,2,1]), [1,2,4,3,5,6])
ap =AllParts(5)
collect(ap)
