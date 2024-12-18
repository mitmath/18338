using SpecialFunctions, Roots, TaylorSeries, LambertW, ProgressLogging

"""
Riemann Siegel theta function
"""
function θ(t)
	return imag.(loggamma(1/4 + t*im/2)) - t*log(π)/2
end

function ψ(p)
	cos(2π*(p^2 - p - 1/16)) / cos(2π*p)
end

taylor_order = 20
N_intervals = 8
a = [(i - 1/2)/N_intervals for i in 1:N_intervals] #center(s) of taylor expansion
t=Taylor1(Float64, taylor_order)
coeffs = [ψ(t+a_i) for a_i in a]

function get_deriv_ψ(k)
    function λ(x)
        s = Int(floor(x*N_intervals) + 1) #interval index
        if s > N_intervals
            s = N_intervals
        end
        
        sum = 0
        for j in range(start=k,stop=taylor_order)
            sum += coeffs[s][j]*(x-a[s])^(j-k)*factorial(j)/factorial(j-k)
        end

        return sum
    end
    return λ
end

dψ = [
    get_deriv_ψ(k) for k in 1:12
]

c_0 = get_deriv_ψ(0)
c_1(p) = -dψ[3](p)/(96π^2)
c_2(p) = dψ[2](p)/(64π^2) + dψ[6](p)/(18432π^4)
c_3(p) = -dψ[1](p)/(64π^2) - dψ[5](p)/(3840π^4) - dψ[9](p)/(5308416π^6)
#c_4(p) = ψ(p)/(128π^2) + 19dψ[4](p)/(24576π^4) + 11dψ[8](p)/(5898240π^6) + dψ[12](p)/(2038431744π^8)

c = Dict(
	0 => c_0,
	1 => c_1,
	2 => c_2,
	3 => c_3,
	#4 => c_4
)

"""
Remainder function in Riemann Siegel formula
"""
function R(t)
		ν = floor(√(t/(2π)))
		μ = √(t/(2π))
		δ = μ - ν
		S = sum(c[k](δ)*μ^(-k) for k in 0:3)
		return (-1)^(ν-1) * μ^-0.5 * S
	end
"""
Riemann Siegel formula: https://mathworld.wolfram.com/Riemann-SiegelFormula.html
"""
function Z(t)
	ν = Int(floor(√(t/(2π))))
	z = 0
	for k in 1:ν
		z += 1/√k * cos(θ(t) - t*log(k))
	end
	return 2*z + R(t)
end

"""
approximation of the n-th Gram point:
https://mathworld.wolfram.com/GramPoint.html
"""
function g(n::Int)
	return 2π * exp(1 + lambertw((8n + 1)/(8ℯ)))
end

"""
computes n-th gram block (g(n), g(n + k)) \\
g(n) must be a "good" gram point
"""
function nth_gram_block(n)
	if (-1)^n * Z(g(n)) <= 0
		error("g(n) is not a good gram point for n = $n")
	end

	k = 1
	while (-1)^(n + k) * Z(g(n + k)) <= 0
		k += 1
	end

	return (n, n+k)
end

"""
input block is a tuple (L, U) \\
finds zeros in the Gram block [g(L), g(U)]
"""
function zeros_in_block(block; tol = 1e-4)

	L, U = block
	if U - L == 1
		return find_zero(Z, (g(L),g(U)), A42(), atol = tol)
	end
	# in the case U - L >= 2, split up the interval and look for sign changes
	g_L = g(L)
	g_U = g(U)
	Δ = g_U - g_L
	ϵ = 4 #size of search grid

	function sign_changes()
		N_splits = (U - L) * ϵ
		search_grid = [g_L + i*Δ/N_splits for i in 0:N_splits]
		Z_vals = [Z(t) for t in search_grid]
		changes = [
			Z_vals[i] < 0 < Z_vals[i + 1] || Z_vals[i] > 0 > Z_vals[i + 1]
			for i in 1:N_splits
		]
		return search_grid, findall(changes)
	end

	grid, idxs = sign_changes()
	while length(idxs) < U - L && ϵ < 500
		ϵ = ϵ*2
		grid, idxs = sign_changes()
	end

	return_arr = Float64[]
	for i in idxs
		ξ = find_zero(Z, (grid[i], grid[i+1]), A42(), atol = tol)
		append!(return_arr, ξ)
		end
	return return_arr
end