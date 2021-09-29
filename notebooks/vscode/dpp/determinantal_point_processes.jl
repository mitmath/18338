### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1b58d214-210a-11ec-091e-f9ac59377f7e
using LinearAlgebra, Combinatorics, Printf, PlutoUI

# ╔═╡ ad9b1c90-3a07-4363-9525-f64e65aa6ed5
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 22907c76-da07-42fe-b641-e1a2977bb8f3
md"""
# Power sets in Julia
"""

# ╔═╡ af37778e-c832-4014-b288-c8e5df86ed5f
md"""
## Slider for N
"""

# ╔═╡ 2c8b8945-a6d9-48c6-a8e4-ada302a00c12
@bind N Slider(1:5, default=4, show_value=true)

# ╔═╡ 27ea925d-336d-4705-87d6-6f82b92feea1
begin
	𝒫 = powerset(1:N)
	collect(𝒫)
end

# ╔═╡ 60b47379-c1e1-4a08-8406-d9c0adf1e59e
md"""
# Det(A+B) Formula
"""

# ╔═╡ c42baf0f-fd04-4fc8-9058-ed3c55e048ab
md"""
``\det(A+B) = \sum_{\mathcal{I}
\in \mathcal{P}} \det(C_\mathcal{I}),``

where ``\mathcal{P}`` is the set of subsets of ``\{1,\ldots,n\}``
and  the mixture of $A$ and $B$ matrix 
``C_\mathcal{I}`` is defined by
$C_\mathcal{I}[:,i] = \left\{ \begin{array}{cc}
 A[:,i] & \text{if } i \notin \mathcal{I} \\
 B[:,i] & \text{if } i \notin \mathcal{I}
\end{array} \right.$
"""

# ╔═╡ 9ae7657d-9d19-48cd-8a66-69be0d1428e8
function mixAandB(A,B,𝓘)
    C = copy(A)
    C[:,𝓘] .= B[:,𝓘]  
    return C
end

# ╔═╡ 0461e492-f7ef-46b2-8d38-d55f25c4bc97
md"""
## Verification:
"""

# ╔═╡ c67b7269-5e43-480c-a7f7-c7bb5bd1bbf7
let  
	### Verification
	A = randn(N,N)
	B = randn(N,N)
	det(A+B), sum(det(mixAandB(A,B,𝓘)) for 𝓘∈𝒫 ) #check identity
end

# ╔═╡ ad82b282-1407-4a24-947b-12789c77707a
md"""
## All mixture matrices for zeros and ones matrix:
"""

# ╔═╡ fd312906-a7e9-48c6-9427-8c4532fc7201
begin
	A = fill(0,N,N)
	B = fill(1,N,N)
	[mixAandB(A,B,𝓘) for 𝓘∈𝒫]
end

# ╔═╡ 6c192d1d-f45a-4ec6-9596-6ad89119e7d1
md"""
# DPP (L definition)
``Pr(𝓘) = \det(L[𝓘,𝓘]) / \det(L+I)``
"""

# ╔═╡ 73725fb6-51fc-45c2-a9c5-0448d6677526
begin
	Y = randn(N,N)
	L = Y'Y
	Pr(𝓘) = det(L[𝓘,𝓘]) / det(L+I)
end

# ╔═╡ 0a5ab519-6f65-4dbb-999a-ee8b83b0141d
md"""
## Explicit probabilities summing to 1
"""

# ╔═╡ 6d9781a8-0d7f-4fdb-a74e-fa5be2aff58a
with_terminal() do 
	for 𝓘 ∈ 𝒫
		@printf("  %10s : %10s \n", 𝓘, Pr(𝓘)) 
    end
	println("--------------------------------------")
	@printf("  %10s : %10s", "sum ",sum(  Pr(𝓘)   for 𝓘 ∈ 𝒫 ) ) 
end

# ╔═╡ 779cdee3-ced4-4da4-9282-2ae855adb40c
md"""
# DPP (K definition)
"""

# ╔═╡ 261be865-1a91-46b3-ab3b-5fc1bc80d231
md"""
$F(\mathcal{I})=Pr(J \supseteq \mathcal{I}) = K \binom{\mathcal{I}}{\mathcal{I}},$
where $K=L(I+L)^{-1}$ or equivalently
$L = K(I-K)^{-1}$
when $I-K$ is invertible.
"""

# ╔═╡ f2ed2cf2-7705-40d3-9959-8f976d1760f7
begin
	K = L/(I+L)
	# To compute Pr(J⊇ 𝓘) just sum the probabilities 
	F( 𝓘 ) = sum( Pr(J)  for J in 𝒫 if J ⊇ 𝓘   )
end

# ╔═╡ 2e4cafc8-c6b7-4569-8d09-1d8e50d469be
with_terminal() do 
	for 𝓘 ∈ 𝒫
		@printf("  %10s : ", 𝓘) 
		@printf("  %10s ", det(K[𝓘,𝓘]))
		@printf("  %10s \n", F( 𝓘 ))
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
Combinatorics = "~1.0.2"
PlutoUI = "~0.7.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "9d8c00ef7a8d110787ff6f170579846f776133a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.4"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "f35ae11e070dbf123d5a6f54cbda45818d765ad2"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.12"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═1b58d214-210a-11ec-091e-f9ac59377f7e
# ╠═ad9b1c90-3a07-4363-9525-f64e65aa6ed5
# ╟─22907c76-da07-42fe-b641-e1a2977bb8f3
# ╟─af37778e-c832-4014-b288-c8e5df86ed5f
# ╠═2c8b8945-a6d9-48c6-a8e4-ada302a00c12
# ╠═27ea925d-336d-4705-87d6-6f82b92feea1
# ╟─60b47379-c1e1-4a08-8406-d9c0adf1e59e
# ╟─c42baf0f-fd04-4fc8-9058-ed3c55e048ab
# ╠═9ae7657d-9d19-48cd-8a66-69be0d1428e8
# ╟─0461e492-f7ef-46b2-8d38-d55f25c4bc97
# ╠═c67b7269-5e43-480c-a7f7-c7bb5bd1bbf7
# ╟─ad82b282-1407-4a24-947b-12789c77707a
# ╠═fd312906-a7e9-48c6-9427-8c4532fc7201
# ╟─6c192d1d-f45a-4ec6-9596-6ad89119e7d1
# ╠═73725fb6-51fc-45c2-a9c5-0448d6677526
# ╟─0a5ab519-6f65-4dbb-999a-ee8b83b0141d
# ╠═6d9781a8-0d7f-4fdb-a74e-fa5be2aff58a
# ╟─779cdee3-ced4-4da4-9282-2ae855adb40c
# ╟─261be865-1a91-46b3-ab3b-5fc1bc80d231
# ╠═f2ed2cf2-7705-40d3-9959-8f976d1760f7
# ╠═2e4cafc8-c6b7-4569-8d09-1d8e50d469be
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
