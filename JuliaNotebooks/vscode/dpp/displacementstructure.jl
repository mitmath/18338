    using LinearAlgebra

    n = 3
    X = rand(n)
    Y = X
    A = rand(n)
    B = rand(n)

    K = (A .* B' .- B.*A') ./ (X.-Y' )

    for i=1:length(A)
        K[i,i] =0.0
    end
    R = K/(I-K)
    Q = (I-K)\A 
    P = (I-K)\B

    RR = (Q*P'-P*Q') ./ (X.-Y'  )
    display(RR)
    display(R)