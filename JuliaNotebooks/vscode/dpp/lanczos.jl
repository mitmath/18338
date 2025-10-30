#Pkg.checkout(“ApproxFun”,”development”)  # need the development branch
using ApproxFun
function lanczos(w,N)
    x = ApproxFun.identity_fun(space(w))

    f1=Fun(1./sqrxt(sum(w)),space(x))

    P = Array(Fun,N + 1)
    β = Array(eltype(w),N)
    γ = Array(eltype(w),N)

    P[1] = f1

    v = x.*P[1]
    β[1] = sum(w.*v.*P[1])

    v = v - β[1]*P[1]
    γ[1] = sqrt(sum(w.*v.^2))

    P[2] = v/γ[1]

    for k = 2:N
        v = x.*P[k] - γ[k-1]*P[k-1]
        β[k] = sum(w.*v.*P[k])
        v = v - β[k]*P[k]
        γ[k] = sqrt(sum(w.*v.^2))
        P[k+1] = v/γ[k]
    end

    P,β,γ
end

w=Fun([1.],GaussWeight())   # exp(-x^2)
n=30
P, β,γ=lanczos(w,n)
norm(γ-sqrt((1:n)/2))  #1.59E-15