### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1b58d214-210a-11ec-091e-f9ac59377f7e
using LinearAlgebra, Combinatorics, Printf, PlutoUI, NamedArrays

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
 B[:,i] & \text{if } i \in \mathcal{I}
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

# ╔═╡ d52d1efe-7248-49e5-aa69-8ea4bdce4a05
md"""
# Möbius Inversion (subsets)
the ζ matrix indicates if A ⊆ B

the μ matrix is the inverse ("inclusion-exclusion principle")
"""

# ╔═╡ a3ec8e0c-e021-4040-9f70-242686866d6d
begin
	ζ = Int.([ A ⊆ B for A∈collect(𝒫), B∈collect(𝒫) ])
	μ =   Int.(inv(ζ))
	labels = string.(𝒫);labels[1] = "∅"
end

# ╔═╡ 967cae80-f548-4c0d-9b01-399f5947e444
ζ

# ╔═╡ 7675b448-9e5e-40f6-afca-43fb66f33ff4
μ 

# ╔═╡ 4cbc0dd1-2d16-47f8-9fbc-ae43e255ca02
NamedArray( ζ , (labels,labels) )

# ╔═╡ 3c870eb7-881d-4636-bdbd-f879407ea854
NamedArray( μ , (labels,labels) )

# ╔═╡ 4fc92447-0a66-44d8-993e-0cb5fceacde8
md"""
## Linking the L and K definitions through the Möbius matrices
"""

# ╔═╡ 3b3aae07-a166-40f9-9973-17e2945c038b
begin
	SetPDF = [det(L[𝓘,𝓘])/det(L+I) for 𝓘∈𝒫] # vector of L definition
	SetCDF = [det(K[𝓘,𝓘])  for 𝓘∈𝒫]  # vector of K definition
	(SetCDF  ≈ ζ  * SetPDF)  && (SetPDF  ≈ μ * SetCDF)
end

# ╔═╡ 70af1efc-4430-41fe-9064-80950d0a1fd1
md"""
# Projection DPP's
What happens to L when K is nearly a projection matrix?
"""

# ╔═╡ 8b673853-cd97-4005-b31f-424c71876bd3
let
	N = 5
	n = 2
	Y = Matrix(qr(randn(N,n)).Q)
    K = Y*Y' .+ .0001*randn(N,N)
	L = K/(I-K)
	
	𝒫 = powerset(1:N)
	
	ζ = Int.([ A ⊆ B for A∈collect(𝒫), B∈collect(𝒫) ])
    SETCDF = [det(K[𝓘,𝓘])  for 𝓘∈𝒫]
	
	with_terminal() do 
	for 𝓘 ∈ 𝒫
		@printf("  %10s : %10s \n", 𝓘, round(det(L[𝓘,𝓘])/det(L+I), digits=3)) 
    end
		
	println("--------------------------------------")
	@printf("  %10s : %10s", "sum  ",sum( det(L[𝓘,𝓘])/det(L+I)   for 𝓘 ∈ 𝒫 ) ) 
		
	
	
	for i in 1:N^2
 	  println( round( (ζ \ SETCDF)[i], digits=3) )
    end
		
		
    end
	

end

# ╔═╡ e7d077dd-c137-4e42-a94e-2e91dc200710
md"""
## Janossy
"""

# ╔═╡ 761439ab-c397-4720-87c9-fb84fe6d1db5
let
	# not fully finished
	N = 8
	n = 2
	Y = Matrix(qr(randn(N,n)).Q)
	
 	K = Y*Y'
	L = zero.(K)
	𝓘 =[2,4,6,8] # Restrict to an "I" 
	
	L[𝓘,𝓘] = K[𝓘,𝓘]
 	L = L/(I-L) # Formula (6) on p.3 
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NamedArrays = "86f7a689-2022-50b4-a561-43c23ac3c673"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
Combinatorics = "~1.0.2"
NamedArrays = "~0.9.6"
PlutoUI = "~0.7.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

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

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "2fd5787125d1a93fbe30961bd841707b8a80d75b"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.9.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "9d8c00ef7a8d110787ff6f170579846f776133a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.4"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "f35ae11e070dbf123d5a6f54cbda45818d765ad2"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.12"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═1b58d214-210a-11ec-091e-f9ac59377f7e
# ╟─ad9b1c90-3a07-4363-9525-f64e65aa6ed5
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
# ╟─d52d1efe-7248-49e5-aa69-8ea4bdce4a05
# ╠═a3ec8e0c-e021-4040-9f70-242686866d6d
# ╠═967cae80-f548-4c0d-9b01-399f5947e444
# ╠═7675b448-9e5e-40f6-afca-43fb66f33ff4
# ╠═4cbc0dd1-2d16-47f8-9fbc-ae43e255ca02
# ╠═3c870eb7-881d-4636-bdbd-f879407ea854
# ╟─4fc92447-0a66-44d8-993e-0cb5fceacde8
# ╠═3b3aae07-a166-40f9-9973-17e2945c038b
# ╟─70af1efc-4430-41fe-9064-80950d0a1fd1
# ╠═8b673853-cd97-4005-b31f-424c71876bd3
# ╟─e7d077dd-c137-4e42-a94e-2e91dc200710
# ╠═761439ab-c397-4720-87c9-fb84fe6d1db5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
