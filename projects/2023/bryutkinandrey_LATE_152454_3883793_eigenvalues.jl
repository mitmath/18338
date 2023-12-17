### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6a9c9a70-aad0-4e13-b90d-76ec85eae514
begin
import Pkg
Pkg.add("Plots")
Pkg.add("DifferentialEquations")
Pkg.add("SpecialFunctions")
Pkg.add("FastGaussQuadrature")
Pkg.add("LinearAlgebra")
Pkg.add("ForwardDiff")
Pkg.add("Interact")
Pkg.add("PlutoUI")
Pkg.add("SpecialPolynomials")
Pkg.add("Polynomials")
Pkg.add("Distributions")
Pkg.add("Combinatorics")
end

# ╔═╡ 0adf0cc6-93ce-11ee-11a9-5b8d798c0eea
using Plots, DifferentialEquations, SpecialFunctions, FastGaussQuadrature, LinearAlgebra, ForwardDiff, Interact, PlutoUI, SpecialFunctions, SpecialPolynomials, Polynomials, Distributions, Combinatorics

# ╔═╡ db918cf9-720b-4c0b-ad61-931225a9e33a
function spacings()
	function PVODE!(du,u,p,t)
		du .= [u[2], -2/t * sqrt(Complex((u[1]-t*u[2])*(t*u[2]-u[1]+u[2]^2))), u[1]/t]
	end

	t0,tn = 5.0,0
	u0 = [-t0/π-(t0/π)^2, -1/π-2*(t0/π), -t0/π-(t0/π)^2/2]
	prob = ODEProblem(PVODE!, u0, (t0,tn))
	sol = solve(prob,Tsit5(), reltol=1e-14, abstol=1e-14)

	print(sol)
	
	#E(t,I) = (t, exp(I))
	
	#PVO(t) = f(t, sol(t)[3]...)[2]

	#TracyWidomPDF_via_Fredholm_Det(s) = ForwardDiff.derivative( t->det(I-K(t)),s)
	
end

# ╔═╡ 821fb6f9-e899-413c-bfa2-498b3677f8c4
begin
	function tracywidomODE!(du,u,p,t) # p is not used
    	du .= [u[2], t*u[1] + 2*u[1]^3, u[4], u[1]^2]
	end
	
	t0,tn = 5.0,-8.0
	u0 = [airy(t0), airy(1, t0), 0, airy(t0)^2]
	prob = ODEProblem(tracywidomODE!, u0, (t0,tn))
	sol = solve(prob,Tsit5(), reltol=1e-14, abstol=1e-14)
	f(t,I,I′) = (t, -I′*exp(-I))
	
	TracyWidomPDF_via_ODE(t) = f(t, sol(t)[[3,4]]...)[2]
end

# ╔═╡ 98e5f1ba-7f15-4bd1-ae8b-6a919d89d493
function DPPSampler()

	function randprojDPP(Y)
    n = size(Y, 2)
    J = fill(0, n)
    for k in 1:n
        p = mean(abs.(Y).^2, dims=2)
        J[k] = rand(Categorical(p[:]))
        Y = ( Y * qr(Y[J[k], :]).Q )[:, 2:end]
    end
    	return sort(J)
	end
	
	function randDPP(Λ, Y)
	   mask = rand.(Bernoulli.(Λ))
	   return randprojDPP(Y[:, mask])#, (mask)
	end

	Γ = gamma
	Hj(j, x) = basis(Hermite, j)(x) # Hermite Polynomial
	ϕ(j, x) = exp(-x^2/2)*Hj(j, x)/( π^(1/4)*sqrt(Γ(j+1))*2^(j/2) )

	Kernel0_Hermite(j, x) = j*ϕ(j,x)^2-sqrt(j*(j+1))*ϕ(j-1,x)*ϕ(j+1,x)
	Kernel_Hermite(j, x, y) = x==y ? Kernel_Hermite(j, x) : sqrt(j/2) * (ϕ(j,x)*ϕ(j-1,y) - ϕ(j-1,x) * ϕ(j,y))/(x-y)
	
	N_DPP = 15
	dx = 0.05
	x = -6:dx:6
	# For example
	K = [Kernel_Hermite(N_DPP, xi, xj) for xi in x, xj in x] * dx

	S, V = eigen(K);
	S[abs.(S).<1e-10] .= 0.0
	S[abs.(S).>1.0] .= 1.0
	r_DPP = [];
	for i in 1:10000
	    append!(r_DPP, randDPP(S, V))
	end

	return histogram(r_DPP, normalized=true, bins=80)
end

# ╔═╡ 0599f826-373a-4ce0-9e2f-041b890303fa
dpp_samples = DPPSampler()

# ╔═╡ 2296f411-4d6b-4296-9877-a44622ab51b8
function GL(n; a=0, b=1) 
    # m points Gauss-Legendre quadrature for [a, b]
    xi, wi = gausslegendre(n)
    return map( t->((a+b)/2 + t*(b-a)/2), xi ), (b-a) * wi / 2
end

# ╔═╡ cdc596e5-5dcd-422a-85a8-bacb1b6176b7
begin
	function K_quadrature(K, s; m=50) 
	    # Creates K matrix according to m quadrature points
	    # Method proposed in (Bornemann, 2010)
	    x, w = GL(m; a=0, b=s)
	    w_sqrt = sqrt.(w)
	    return (w_sqrt*w_sqrt') .* map(t->K(t...), [(xi,xj) for xi in x, xj in x]) 
	end

	function K_quadrature_bessel(K, α, s; m=50) 
	    # Creates K matrix according to m quadrature points
	    # Method proposed in (Bornemann, 2010)
	    x, w = GL(m; a=0, b=last(s))
	    w_sqrt = sqrt.(w)
	    return (w_sqrt*w_sqrt') .* map(t->K(α,s,t...), [(xi,xj) for xi in x, xj in x]) 
	end
end

# ╔═╡ eb143922-4a7e-462b-a2c6-747eaaaf3ac8
begin
K_sinc(x, y)    = x==y ? 1.0 : sin( π * (x - y) ) / ( π * (x - y) )
K0_sinc(x) = sin( π*x )/ (π*x)

F₀(s; m=30) = det(I - K_quadrature(Kp_sinc, s; m=m))
F = plot(xs, F₀.(xs), xlabel="s", leg=false, lw=3, color=:black, ylims=(0,1), xlims=(0,4), size=(400, 300))
plot(F)
end

# ╔═╡ 58159a73-51fd-4619-9ce0-a95370879464
begin
"""
HEpdf  : pdf of the smallest eigenvalue, hard-edge
HE2pdf : pdf of the second smallest eigenvalue, hard-edge
"""
	Ja = besselj
	K_Bess(α, x, y) = x==y ? ( Ja(α, √x)^2 - Ja(α+1, √x) * Ja(α-1, √x) ) / 4 :
	                                 ( Ja(α+1, √x) * √x * Ja(α, √y) - √y * Ja(α+1, √y) * Ja(α, √x) ) / ( 2 * (x - y) )
	# Again, everything what has a p in the name is conditional and can be removed
	
	function HE2pdf(α, s; m=50)
	    K = K_quadrature_bessel(K_Bessp, α, s; m=m)
	    L = (I-K)\K
	    return K_Bess(α,s,s)*tr(L)/det(I+L)
	end
end

# ╔═╡ 8561c70c-1faa-4dc1-99bf-aeecc0519e65
begin
	airy_kernel(x, y) = x==y ? (airyaiprime(x))^2 - x * (airyai(x))^2 :
	           (airyai(x) * airyaiprime(y) - airyai(y) * airyaiprime(x)) / (x - y)
	
	ϕ(ξ, s) =  s + 10*tan(π*(ξ+1)/4) # Transformation from [-1,1] to (s,∞)
	ϕ′(ξ) = (5π/2)*(sec(π*(ξ+1)/4))^2
	K(ξ,η,s) = sqrt(ϕ′(ξ) * ϕ′(η)) * airy_kernel(ϕ(ξ,s), ϕ(η,s))
	
	function K(s , n=100) 
	    nodes,weights = gausslegendre(n)
	    Symmetric( K.(nodes',nodes,s) .* (√).(weights) .* (√).(weights'))
	end
	
	TracyWidomPDF_via_Fredholm_Det(s) = ForwardDiff.derivative( t->det(I-K(t)),s)
end

# ╔═╡ 6c9adf92-12e0-4119-a269-5de79994586b
function conditionalDPP_pdf()
	
	function K_quadrature(K, s; m=50) 
	    x, w = GL(m)
	    w_sqrt = sqrt.(w)
	    return (w_sqrt*w_sqrt') .* map(t->K(s, t...), [(xi,xj) for xi in x, xj in x]) 
	end
	
	Ai, Aip = airyai, airyaiprime
	ϕ(x,s)    = s+10*tan(π*x/2) # See (Bornemann, 2010) pdf page 29
	ϕp(x)     = 5π*sec(π*x/2)^2
	K_Ai(x,y)      = x==y ? (Aip(x))^2-x*(Ai(x))^2 : (Ai(x)*Aip(y)-Ai(y)*Aip(x))/(x-y)
	
end

# ╔═╡ 6d68ec96-fb7a-4600-a543-2f1a1339b994
conditionalDPP_pdf()

# ╔═╡ 8b1ecc85-12f9-480b-82e8-76e2303b2a42


# ╔═╡ 6fe64efe-06e2-471f-9f0f-c93c7525a726
begin
	t = 10_000
	e = fill(0.0,t)
	n = 100
	@time for i = 1:t
	    A = ( randn(n,n) + im * randn(n,n))
	    e[i] = eigmax(Hermitian( A+A'))   ## Random Matrices from the GUE, and take eigmax
	end
	e .=   (e/2 .- 2*√n) * n^(1/6) # normalize
	h = histogram(e, normalized=true, bins=80)
end

# ╔═╡ b201d03b-ebb4-46c2-ae16-6359044579ee
begin
	pODE = plot(sol, vars=(f,0,3,4), flip=false, legend=false)
	pDet = plot( -8:.1:5, TracyWidomPDF_via_Fredholm_Det)
	pHist = plot(h,color=:yellow, label="Experiment")
	pDPP = plot(dpp_samples, color=:red, label="Experiment")
	plot(pODE, pDet, pHist, pDPP)
end

# ╔═╡ af2c4da8-f735-4b84-ab24-48ead8b5d16e
begin
plot(h,color=:yellow, label="Experiment")
plot!(-8:.1:5, TracyWidomPDF_via_ODE, lw=3, color=:red, label="PDF")
plot!(-8:.1:5, TracyWidomPDF_via_Fredholm_Det, label="Fredholm" )
end

# ╔═╡ 6125b43d-9157-43b8-8aa7-e273cfc8c3dc
function randDPP()
	function randDPP_P(Kernel)
	   K = copy(Kernel)
	   N = size(K,2) 
	   𝓘 =  Int[]
	   for i=1:N
	     rand(Bernoulli(K[i,i])) ?  push!(𝓘 ,i) :  K[i,i] -= 1
	     K[i+1:N,i+1:N] -= K[i,i+1:N]*K[i+1:N,i]'/K[i,i]
	   end
	   return 𝓘 
	end
	
	N = 3

	L = randn(N,N)
	L *= L'
	K_main = L/(L+I)
	
	
	hist = Dict( 𝓘=>0 for 𝓘 ∈ powerset(1:N) )
	t = 10000
	for i=1:t
	    hist[randDPP_P(K_main)] += 1
	end
	
	println("$t trials N=$N")
	println("Expmnt Theory")
	
	pp(𝓘) = length(𝓘)>0 ? 𝓘 : "∅" # pretty print the empty set
	
	for 𝓘∈powerset(1:N)
	   println(round( hist[𝓘]/t,digits=3)," ",round.(det(L[𝓘,𝓘])/det(L+I),digits=3), " ",pp(𝓘))
	end
end

# ╔═╡ 20b2099a-b570-4e6c-b0ac-33402d1f35d1
randDPP()

# ╔═╡ d057d97e-6123-428f-8041-1db4057c941f
function generate_beta_tridiagonal_matrix(n::Int, β::Float64)
    # Diagonal elements: N(0, 2)
    main_diag = [rand(Normal(0, sqrt(2))) for _ in 1:n]

    # Off-diagonal elements: Chi distribution
    off_diag = [rand(Chisq((n-i)*β)) for i in 1:n-1]

    # Creating the Tridiagonal matrix
    H = Tridiagonal(off_diag, main_diag, off_diag)
    return H
end

# ╔═╡ 656fff8a-d925-4588-9f03-081fbaa5e217
function biggest_eigenvalue_gaussian(β)
	t = 100
	e = fill(0.0,t)
	n = 100
	eigenvalues = []
	eigenvalues_normal = []
	A = randn(n,n)
	for i = 1:t
	    e[i] = eigmax(generate_beta_tridiagonal_matrix(n,β))
		eigenvals = eigen(generate_beta_tridiagonal_matrix(n, β)).values
        append!(eigenvalues, eigenvals)

		
		eigenvals_normal = eigvals(Symmetric( A+A')/2)
        append!(eigenvalues, eigenvals)
		append!(eigenvalues_normal, eigenvals_normal)
	end
	
	h_eig = histogram(eigenvalues, normalized=true, bins=200, legend=false, xlabel="Eigenvalue", ylabel="Frequency", title="Histogram of Eigenvalues")
	h_eig_normal = histogram(eigenvalues_normal, normalized=true, bins=200, legend=false, xlabel="Eigenvalue", ylabel="Frequency", title="Histogram of Eigenvalues")
	e .=   (e/2 .- 2*√n) * n^(1/6) # normalize
	h = histogram(e, normalized=true, bins=200)
	plot(h, h_eig,h_eig_normal)
end

# ╔═╡ bfdffd52-b052-49e0-a0ab-d4dc502f2d0a
 biggest_eigenvalue_gaussian(2.0)

# ╔═╡ 84b8ee2a-7927-44b1-b8cb-6b8ba850374f
begin
	function generate_B_bidiagonal_matrix(n::Int, β::Float64, a::Float64 = ceil(β/2 * (n+1)))
	    # Ensuring 'a' satisfies the condition
	    a = max(a, ceil(β/2 * (n+1)))
	
	    # Diagonal elements: chi-distributed with degrees of freedom 2a - β(i-1)
	    main_diag = [rand(Chisq(2 * a - β * (i - 1))) for i in 1:n]
	
	    # Sub-diagonal elements: chi-distributed with degrees of freedom β(n-i)
	    sub_diag = [rand(Chisq(β * (n - i))) for i in 1:n-1]
	
	    # Creating the Bidiagonal matrix
	    B = Bidiagonal(main_diag, sub_diag, :L)
	    return B
	end
	
	function generate_L_bidiagonal_matrix(n::Int, β::Float64, a::Float64 = ceil(β/2 * (n+1)))
	    B = generate_B_bidiagonal_matrix(n, β, a)
	    L = Tridiagonal(B * B')
	    return L
	end
	
	# Example usage
	#n = 5
	#β = 1.0
	#generate_L_bidiagonal_matrix(n, β)  # 'a' is optional and calculated if not provided

end

# ╔═╡ 11870736-95d1-425e-9ee5-dc1115abe8aa
function biggest_eigenvalue_laguerre(β)
	t = 100
	e = fill(0.0,t)
	eigenvalues = []
	n = 100
	for i = 1:t
	    e[i] = eigmax(generate_L_bidiagonal_matrix(n,β))
		eigenvals = eigen(generate_L_bidiagonal_matrix(n, β)).values
        append!(eigenvalues, eigenvals)
	end
	#e .=   (e/2 .- 2*√n) * n^(1/6) # normalize
	h_eig = histogram(eigenvalues, normalized=true, bins=200, legend=false, xlabel="Eigenvalue", ylabel="Frequency", title="Histogram of Eigenvalues")
	h = histogram(e, normalized=true, bins=200)
	plot(h,h_eig)
end

# ╔═╡ fc9b6ab9-92bc-452c-9767-dc90dd926945
biggest_eigenvalue_laguerre(2.0)

# ╔═╡ 1361001d-be51-4f39-b2c4-0e9bcb90d5cb
begin
function generate_J_matrix(n::Int, β::Float64, a::Float64, b::Float64)
    # Generating beta-distributed values and calculating angles
	    c = [sqrt(rand(Beta(β/2*(a+i), β/2*(b+i)))) for i in 1:n]
	    c_prime = [sqrt(rand(Beta(β/2*i, β/2*(a+b+1+i)))) for i in 1:n-1]
	    s = sqrt.(c)
	    s_prime = sqrt.(c_prime)
	
	    # Constructing the sub-matrices as Bidiagonal
	    B11 = Bidiagonal(reverse([c[1:end-1] .* s_prime; c[end]]), -reverse(s[2:end] .* c_prime), :U)
	    B12 = Bidiagonal(reverse([s[1];s[2:end] .* s_prime]), reverse(c[1:end-1] .* c_prime), :L)
	    B21 = Bidiagonal(-reverse([s[1:end-1] .* s_prime; s[end]]), -reverse(c[2:end] .* c_prime), :U)
	    B22 = Bidiagonal(reverse([c[1]; c[2:end] .* s_prime]),-reverse(s[1:end-1] .* c_prime), :L)
	
	    # Combining sub-matrices into the final matrix
	    upper_half = [B11 B12]
	    lower_half = [B21 B22]
	    J = [upper_half; lower_half]
	    return J
	end
	
	# Example usage
	#N = 5
	#β = 1.0
	#a = 0.5
	#b = 0.5
	#generate_J_matrix(N, β, a, b)
end

# ╔═╡ Cell order:
# ╠═6a9c9a70-aad0-4e13-b90d-76ec85eae514
# ╠═0adf0cc6-93ce-11ee-11a9-5b8d798c0eea
# ╠═db918cf9-720b-4c0b-ad61-931225a9e33a
# ╠═821fb6f9-e899-413c-bfa2-498b3677f8c4
# ╠═98e5f1ba-7f15-4bd1-ae8b-6a919d89d493
# ╠═0599f826-373a-4ce0-9e2f-041b890303fa
# ╠═2296f411-4d6b-4296-9877-a44622ab51b8
# ╠═cdc596e5-5dcd-422a-85a8-bacb1b6176b7
# ╠═eb143922-4a7e-462b-a2c6-747eaaaf3ac8
# ╠═58159a73-51fd-4619-9ce0-a95370879464
# ╠═8561c70c-1faa-4dc1-99bf-aeecc0519e65
# ╠═6c9adf92-12e0-4119-a269-5de79994586b
# ╠═6d68ec96-fb7a-4600-a543-2f1a1339b994
# ╠═8b1ecc85-12f9-480b-82e8-76e2303b2a42
# ╠═6fe64efe-06e2-471f-9f0f-c93c7525a726
# ╠═b201d03b-ebb4-46c2-ae16-6359044579ee
# ╠═af2c4da8-f735-4b84-ab24-48ead8b5d16e
# ╠═6125b43d-9157-43b8-8aa7-e273cfc8c3dc
# ╠═20b2099a-b570-4e6c-b0ac-33402d1f35d1
# ╠═d057d97e-6123-428f-8041-1db4057c941f
# ╠═656fff8a-d925-4588-9f03-081fbaa5e217
# ╠═bfdffd52-b052-49e0-a0ab-d4dc502f2d0a
# ╠═84b8ee2a-7927-44b1-b8cb-6b8ba850374f
# ╠═11870736-95d1-425e-9ee5-dc1115abe8aa
# ╠═fc9b6ab9-92bc-452c-9767-dc90dd926945
# ╠═1361001d-be51-4f39-b2c4-0e9bcb90d5cb
