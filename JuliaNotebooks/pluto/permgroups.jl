### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ c79f0fda-4fc0-11ed-330a-4d48d5a6b34d
using LinearAlgebra, Combinatorics, Random, SparseArrays

# ╔═╡ 3c5ed858-aea8-4633-80ca-876ca1eb3295
md"# Make a group"

# ╔═╡ 7eaba50e-4afe-4473-9b53-eca602a8914a
begin
  Base.ComposedFunction(π::Vector, σ::Vector)  = π[σ]
  inv(π::Vector) = sortperm(π)
end

# ╔═╡ d1d4eafd-dc3f-4319-9e94-a90facf7600e
n = 3

# ╔═╡ 9f7d94dc-26d4-4757-b601-00849c4d7b50
π = randperm(n)

# ╔═╡ 70c85fc8-926b-4cdf-9464-853252e41381
σ = randperm(n)

# ╔═╡ e4d66bc2-4800-4400-8801-94191151d231
τ = randperm(n)

# ╔═╡ 1a88137e-d009-4e6e-9222-a6991963b0e6
(π ∘ σ) ∘ τ 

# ╔═╡ 6cdd7fc8-748b-44a8-a61f-1258b70d36a9
π ∘ (σ ∘ τ)

# ╔═╡ bdd2b1c5-5b35-4a59-b201-487bc2e1ba91
inv(π) ∘ π

# ╔═╡ 24c1e3b7-a9ba-4022-966f-e279dcaed50f
md"# Permutation Representation"

# ╔═╡ 6292248b-20fe-42fc-a0ce-e3893cfdf834
ρ(π::Vector) = sparse(I(length(π))[:,π])

# ╔═╡ 45b77bd8-7e99-4d6a-aa14-9b7aede01826
ρ(π )

# ╔═╡ 269f5baf-2dce-4248-9380-cc4d2149dfcd
ρ( π ∘ σ ) , ρ(π) * ρ(σ)

# ╔═╡ 5e65b263-1e25-4333-8624-4999d14c0a8d
md"# Rotated Permutations"

# ╔═╡ 6ebed54c-a004-40da-8675-f831c5c148a6
Q = qr(ones(n,1)).Q[:,1:n]

# ╔═╡ a795c678-a89c-467c-ab51-be6bccdbfa24
ρQ( π ) =  Q * ρ(π) * Q

# ╔═╡ f3959f0c-f2de-4bb2-b5a2-34e2e64684da
 [  sparse(round.(ρQ( π ),digits=3)) for π ∈ permutations(1:n) ]

# ╔═╡ dad6007b-a410-4513-97f7-89d62880f8dc
md"# Projected Permutations"


# ╔═╡ e33ea769-eb77-4e4a-9666-9606f9e3b572
begin
  Π = Q[:,2:n]
  ρΠ(π) = Π' * ρ(π) * Π
  [  sparse(round.(ρΠ( π ),digits=3)) for π ∈ permutations(1:n) ]
end

# ╔═╡ fcb12973-ef9a-4408-8900-9fc09ed401e9
tr.(ρΠ( π ) for π ∈ permutations(1:n))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
Combinatorics = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-rc4"
manifest_format = "2.0"
project_hash = "4e386a005d7a5a6f799f5ddf368ff649c3e6ba25"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╠═c79f0fda-4fc0-11ed-330a-4d48d5a6b34d
# ╟─3c5ed858-aea8-4633-80ca-876ca1eb3295
# ╠═7eaba50e-4afe-4473-9b53-eca602a8914a
# ╠═d1d4eafd-dc3f-4319-9e94-a90facf7600e
# ╠═9f7d94dc-26d4-4757-b601-00849c4d7b50
# ╠═70c85fc8-926b-4cdf-9464-853252e41381
# ╠═e4d66bc2-4800-4400-8801-94191151d231
# ╠═1a88137e-d009-4e6e-9222-a6991963b0e6
# ╠═6cdd7fc8-748b-44a8-a61f-1258b70d36a9
# ╠═bdd2b1c5-5b35-4a59-b201-487bc2e1ba91
# ╠═24c1e3b7-a9ba-4022-966f-e279dcaed50f
# ╠═6292248b-20fe-42fc-a0ce-e3893cfdf834
# ╠═45b77bd8-7e99-4d6a-aa14-9b7aede01826
# ╠═269f5baf-2dce-4248-9380-cc4d2149dfcd
# ╠═5e65b263-1e25-4333-8624-4999d14c0a8d
# ╠═6ebed54c-a004-40da-8675-f831c5c148a6
# ╠═a795c678-a89c-467c-ab51-be6bccdbfa24
# ╠═f3959f0c-f2de-4bb2-b5a2-34e2e64684da
# ╠═dad6007b-a410-4513-97f7-89d62880f8dc
# ╠═e33ea769-eb77-4e4a-9666-9606f9e3b572
# ╠═fcb12973-ef9a-4408-8900-9fc09ed401e9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
