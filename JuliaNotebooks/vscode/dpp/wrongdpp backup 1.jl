function wrongDPP(K)
    n = size(K,2)
    𝓘 = fill(0,n)
    for k=1:n
        p = mean(abs.(Y).^2, dims=2)
        𝓘[k] = rand(Categorical(p[:]))
        Y=(Y*qr(Y[𝓘[k],:]).Q )[:,2:end] 
        display(Y[𝓘[k],:]) 
    end
    return(sort(𝓘))
end