### A Pluto.jl notebook ###
# v0.11.10

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

# ╔═╡ 5dd6cbdc-f2bc-11ea-296b-a3f2aac6c1e1
using LinearAlgebra, Plots, PlutoUI

# ╔═╡ 6c267c96-f2bc-11ea-0764-7dc3e141d681
begin
	# semicircle.jl
	# Algorithm 1.2 of Random Eigenvalues by Alan Edelman
	
	# Experiment:  Sample random symmetric Gaussian matrices
	# Plot:        Histogram of the eigenvalues
	# Theory:      Semicircle as n->infinity
	
	## Parameters
	n = 1000   # matrix size
	t = 1      # trials
	
	## Experiment

	
	function GOElike_eigvals(n, rand_func) 
	    A = [rand_func() for i=1:n, j=1:n] # draw n by n matrix 
	    S = Symmetric((A + A') / 2)    # symmetrize matrix
	    return eigvals(S)
	end
	
	function sample(n, t, rand_func)         # run t trials with matrix size n
	    v = Float64[]                  # initialize result vector
	    for i = 1:t
	        append!(v, GOElike_eigvals(n, rand_func))
	    end
	    v ./= sqrt(n / 2)              # scale eigenvalues
	    return v
	end
	
	## Plot
	function semicircle(n, t, rand_func = randn) # default argument
	    v = sample(n, t, rand_func)                        # run experiment
	    histogram(v, bins=100, normed = true)                  # create histogram
	    plot!(-2:0.01:2, x -> √(4 - x^2)/(2π), c=:red, lw=2, aspectratio=6.5, leg=false)   # plot semicircle
	    xlims!(-3, 3)
		 ylims!(0,.5) ;                                  # set limits for x axis in plot
	end
end

# ╔═╡ 7904a030-f2bc-11ea-1a46-3b4588e10ead
md"""
 $(@bind nn Slider(1:100,show_value=true))
 $(@bind tt Slider(1:10000,show_value=true))
"""

# ╔═╡ 83755f36-f2bc-11ea-0080-6b845cdae12c
semicircle(nn, tt, ()->randn())

# ╔═╡ 9d0f2de4-f2bc-11ea-0020-1beccc39360f
semicircle(nn, tt,  ()->(rand()-0.5)*sqrt(12))

# ╔═╡ Cell order:
# ╠═5dd6cbdc-f2bc-11ea-296b-a3f2aac6c1e1
# ╠═6c267c96-f2bc-11ea-0764-7dc3e141d681
# ╠═7904a030-f2bc-11ea-1a46-3b4588e10ead
# ╠═83755f36-f2bc-11ea-0080-6b845cdae12c
# ╠═9d0f2de4-f2bc-11ea-0020-1beccc39360f
