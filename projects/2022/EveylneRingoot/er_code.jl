using LinearAlgebra, Distributions, Random, StaticArrays, SparseArrays, Plots, StatsPlots, DelimitedFiles, Plots.Measures, CUDA, BenchmarkTools
using Combinatorics, DifferentialEquations, SpecialFunctions, FastGaussQuadrature, ForwardDiff, Interact, Statistics, StatsBase


function bidiag(n,m)
    (m>n) && error("n must be larger or equal than m");
    return Float32.([diagm([rand(Chi(k)) for k in n:-1:(n-m+1)]); zeros(1,m)] .+ [zeros(1,m); diagm([rand(Chi(k)) for k in m:-1:1])])
end

function bidiag_trunc(n,m,select)
    (m>n) && error("n must be larger or equal than m");
    (select>m) && error("select must be smaller than m");
    return Float32.(diagm([rand(Chi(k)) for k in n:-1:(n-select+1)]) .+ [zeros(1,select); [diagm([rand(Chi(k)) for k in m:-1:m-select+2]) zeros(select-1,1)]])
end

#largest singular values
n_values=[40, 60, 100, 140, 200, 280, 400]
Truncating_amount=[0.05 ,0.1, 0.2]
M_larger_than_n_value=1

trials=100
myplot=plot(xlabel="m=10n", title="Absolute error of truncated svd estimation \n over "*string(trials)*" trials with m="*string(M_larger_than_n_value)*"n", legend = :bottomleft, size=(400,320), legendfontsize=7, yaxis=:log, titlefontsize=10)

j=1
for  red in Truncating_amount
    acc=zeros(trials,length(n_values))
    i=1
    for n in n_values
        for t in 1:trials
            testmatrix=bidiag(M_larger_than_n_value*n,n)
            testmatrix_trunc=testmatrix[1:Int(red*n),1:Int(red*n)]
            real_svd=svdvals(testmatrix)[1]
            approx_svd=svdvals(testmatrix_trunc)[1]
            acc[t,i]=abs(approx_svd-real_svd)./real_svd
        end
        i=i+1
    end
    myplot=plot!(string.(transpose((n_values.*M_larger_than_n_value)[:,:])), acc,  seriestype = :violin, fillalpha=0.5, linealpha=0.5, markeralpha=0.5, color=palette(:tab10)[j], 
            labels=permutedims(["Portion of matrix used = "*string(red),["" for i=1:length(n_values)-1]...]))
    j=j+1
end


################################################################################################
#-------------------------Sensitivity of accuracy to truncation-----------------------------
################################################################################################

n_values=unique([round(Int, 10^i) for i=0:0.2:3.4])
trials=150
min_n_for_k=[1,6,11]
k_values=[1,10,100]
alpha_values=[1,10,100,1000]

#Im pre-allocating some stuff on the GPU so it doesn't require any allocations during the calculations
seriesize_gpu=15 
N_switch_to_gpu=150
truncated_n_values=sort(unique(round.(Int,n_values ./n_values')))
truncated_n_values_gpu=truncated_n_values[truncated_n_values.>N_switch_to_gpu]
n_values_gpu=n_values[n_values.>N_switch_to_gpu]
inputs_gpu_square=[CUDA.zeros(svdsize,svdsize)  for _ in 1:seriesize_gpu for svdsize in truncated_n_values_gpu];
inputs_gpu_square=reshape(inputs_gpu_square, length(truncated_n_values_gpu), seriesize_gpu);
inputs_gpu_squareplusone=[CUDA.zeros(svdsize+1,svdsize)  for _ in 1:seriesize_gpu for svdsize in n_values_gpu];
inputs_gpu_squareplusone=reshape(inputs_gpu_squareplusone, length(n_values_gpu), seriesize_gpu);
outputs_gpu=[CUDA.zeros(svdsize,1) for _ in 1:seriesize_gpu for svdsize in truncated_n_values_gpu];
outputs_gpu=reshape(outputs_gpu, length(truncated_n_values_gpu), seriesize_gpu);

for (alpha_i, alpha) in enumerate(alpha_values)
    println("Starting alpha "*string(alpha))
    result_avg=[zeros(length(n_values),length(n_values)) for _ in 1:length(k_values)]
    result_max=[zeros(length(n_values),length(n_values)) for _ in 1:length(k_values)]
    for (n_i,n) in enumerate(n_values)
        println("n is equal to "*string(n));
        testmatrices=[bidiag(alpha*n,n) for _ in 1:trials]
        real_svd=compute_svd_series(testmatrices, outputs_gpu, truncated_n_values_gpu, inputs_gpu_squareplusone, n_values_gpu,
                            N_switch_to_gpu,k_values, trials, seriesize_gpu)
        acc=[zeros(trials, length(k_values)) for _ in 1:length(n_values)]
        for (p_inv_i,p_inv) in enumerate(n_values[1:n_i])
            (round(Int, n./p_inv)<1) && continue;
            testmatrices_trunc=[testmatrix[1:round(Int,n./p_inv),1:round(Int, n./p_inv)] 
                                    for testmatrix in testmatrices]
            approx_svd=compute_svd_series(testmatrices_trunc, outputs_gpu, truncated_n_values_gpu, inputs_gpu_square, truncated_n_values_gpu,
                                        N_switch_to_gpu,k_values, trials, seriesize_gpu)
            No_k_considered=size(approx_svd,1)
            acc[end+1-p_inv_i][:,1:No_k_considered]=max.(abs.((real_svd[1:No_k_considered,:]-approx_svd)./real_svd[1:No_k_considered,:])',1e-7)
        end
        total_mean=reduce(vcat,mean.(acc,dims=1))
        total_max=reduce(vcat,maximum.(acc,dims=1))
        for (k_i,k) in enumerate(k_values)
            result_avg[k_i][n_i,:]=total_mean[:,k_i]
            result_max[k_i][n_i,:]=total_max[:,k_i]
        end
    end
    for (k_i,k) in enumerate(k_values)
        writedlm("svd_avg_mn"*string(alpha)*"k"*string(k) * ".txt" , result_avg[k_i])
        writedlm("svd_max_mn"*string(alpha)*"k"*string(k) * ".txt" , result_max[k_i])
    end
end

################################################################################################
#-------------------------Parallel function for comuting svd faster-----------------------------
################################################################################################
function do_svd!(c,a)
    c .= svdvals!(a);
end

function compute_svd_series(inputs, all_outputs_gpu, n_values_output, zeros_gpu, zeros_n_values,n_gpu,k_values, n_inputs, batchsize)
    current_size=minimum(size(first(inputs)))
    current_k_values=k_values[k_values.<=current_size]
    if current_size<n_gpu
        output_total=svdvals.(inputs)
        output=[output[current_k_values] for output in output_total]
    else
        output=[]
        outputs_gpu=all_outputs_gpu[findfirst(x->x==current_size,n_values_output),:]
        inputs_gpu=zeros_gpu[findfirst(x->x==current_size,zeros_n_values),:]
        for beginindex in 1:batchsize:n_inputs
            for (i,input) in enumerate(inputs[beginindex:min(beginindex+batchsize-1,end)])
                inputs_gpu[i][diagind(input)]=(input[diagind(input)])
                inputs_gpu[i][diagind(input,-1)]=(input[diagind(input,-1)])
            end
            @sync for i in eachindex(inputs_gpu)
                Threads.@spawn begin
                    do_svd!(outputs_gpu[i], inputs_gpu[i])
                    CUDA.synchronize()
                end
            end
            append!(output, [Array(output[current_k_values]) for output in outputs_gpu])
        end
    end
    return reduce(hcat,output)
end

function compute_single_svd_series(inputs,inputs_gpu, outputs_gpu, n_gpu,k, n_inputs, batchsize)
        if size(first(inputs))[1]<n_gpu
        output_total=svdvals.(inputs)
        output=[output[k] for output in output_total]
    else
        output=[]
        for beginindex in 1:batchsize:n_inputs
            for (i,input) in enumerate(inputs[beginindex:min(beginindex+batchsize-1,end)])
                inputs_gpu[i][diagind(input)]=(input[diagind(input)])
                inputs_gpu[i][diagind(input,-1)]=(input[diagind(input,-1)])
            end
            @sync for i in eachindex(inputs_gpu)
                Threads.@spawn begin
                    do_svd!(outputs_gpu[i], inputs_gpu[i])
                    CUDA.synchronize()
                end
            end
            append!(output, [output[k] for output in outputs_gpu])
        end
    end
    return output
end

################################################################################################
#-------------------------Heatmaps-----------------------------
################################################################################################


fontsize=14
ticks=[1,10,100,1000]
cticks=[-2,-4,-6]

pyplot()
allplots=[]
for (trialsummarytype_i,trialsummarytype) in enumerate(["avg","max"])
    for (alpha_i, alpha) in enumerate(alpha_values)
        for (k_i,k) in enumerate(k_values)
            n_values_k=n_values[min_n_for_k[k_i]:end]
            xlims=(minimum(n_values_k),maximum(n_values_k))
            temp_ticks=ticks[k_i:end]

            result=readdlm("svd_"*trialsummarytype*"_mn"*string(alpha)*"k"*string(k) * ".txt" )
            result=log10.(result[min_n_for_k[k_i]:end,min_n_for_k[k_i]:end])
            plot(size=(1000,600))
            myplot=plot!(n_values_k, reverse(1 ./n_values_k), result', 
                levels=[6,6,5,5][alpha_i], 
                xaxis=:log10, yaxis=:log10, 
                seriestype=:contour,  c=:viridis,
                fill=false,lw=2,
                xticks=(temp_ticks,string.(temp_ticks).*"\n".*string.(temp_ticks.*alpha)), 
                yticks=(1 ./temp_ticks[end:-1:1],string.(1 ./temp_ticks[end:-1:1])) , 
                xtickfontsize=fontsize,ytickfontsize=fontsize,
                ylabel="p",xlabel="Matrix size", yguidefontrotation=-90, 
                yguidefontsize=fontsize, xguidefontsize=fontsize,
                margin=1cm,
                title="k = "*string(k)*", relative error "*["averaged over ","as maximum of "][trialsummarytype_i]*string(trials)*" trials", 
                titlefontsize=fontsize,
                colorbar_tickfontsize=fontsize-1, 
                colorbar_ticks=(cticks, string.(round.(10. .^cticks,digits=6)))
                )
            myplot=plot!(n_values_k, reverse(1 ./n_values_k), result', 
                xaxis=:log10, yaxis=:log10, 
                seriestype=:heatmap, c=:viridis,
                fillalpha=0.5, 
                legend=:none, 
                xlims=xlims, ylims=reverse(1 ./xlims)
                )
            myplot=annotate!(2000, 1 ./xlims[end],
                text("m\nn",:left,:top,fontsize+1)
                )
            myplot=annotate!(xlims[2],[1.3,0.12, 0.0115][k_i],
                text("Relative error",:left,fontsize-1)
            )
            push!(allplots,myplot)
            savefig("svd_"*trialsummarytype*"_mn"*string(alpha)*"k"*string(k)*".png")
        end
    end
end



################################################################################################
#-------------------------Verification of approximate formula-----------------------------
################################################################################################

#Hint when trying this at home: clear GPU memory first
trials=20
n=10000
p=6.4*(n^(-0.8))
n_select=round(Int,n*p)
acc=zeros(trials)
for trial in 1:trials
    GC.gc(true)
    CUDA.reclaim()
    A=CuArray(bidiag(n,n))
    A_trunc=deepcopy(A[1:n_select,1:n_select])

    svd_A=svdvals!(A)
    svd_A_trunc=svdvals!(A_trunc)
    svd_A1=svd_A[1]
    svd_A_trunc1=svd_A_trunc[1]
    acc[trials]=(svd_A1-svd_A_trunc1)./svd_A1
end
mean(acc)
################################################################################################
#-------------------------limit distribution-----------------------------
################################################################################################

trials=150
n=400
p=6.4*(n^(-0.8))
n_select=round(Int,n*p)
n_select2=round(Int,n*p/2)
n_select3=round(Int,n*p/5)
n_select4=round(Int,n*p/10)
batchsize=10

trial_matrices=[bidiag_trunc(n,n,n) for _ in 1:trials]
inputs_gpu_n=[CUDA.zeros(n,n)  for _ in 1:batchsize];
outputs_gpu_n=[CUDA.zeros(n,1) for _ in 1:batchsize];
svds_n=compute_single_svd_series(trial_matrices,inputs_gpu_n, outputs_gpu_n, 150,1, trials, batchsize)
svds_scaled_n=((Float64.(svds_n)).^2 .-4 *(n))./(2(2(n)).^(1/3))

trial_matrices_trunc=[bidiag_trunc(n,n,n_select) for _ in 1:trials]
inputs_gpu_2=[CUDA.zeros(n_select,n_select)  for _ in 1:batchsize];
outputs_gpu_2=[CUDA.zeros(n_select,1) for _ in 1:batchsize];
svds_2=compute_single_svd_series(trial_matrices_trunc,inputs_gpu_2, outputs_gpu_2, 150,1, trials, batchsize)
svds_scaled_2=((svds_2) .^2 .-4 *(n))./(2(2(n)).^(1/3))

trial_matrices_trunc=[bidiag_trunc(n,n,n_select2) for _ in 1:trials]
inputs_gpu_3=[CUDA.zeros(n_select2,n_select2)  for _ in 1:batchsize];
outputs_gpu_3=[CUDA.zeros(n_select2,1) for _ in 1:batchsize];
svds_3=compute_single_svd_series(trial_matrices_trunc,inputs_gpu_3, outputs_gpu_3, 150,1, trials, batchsize)
svds_scaled_3=((svds_3) .^2 .-4 *(n))./(2(2(n)).^(1/3))

trial_matrices_trunc=[bidiag_trunc(n,n,n_select3) for _ in 1:trials]
inputs_gpu_3=[CUDA.zeros(n_select3,n_select3)  for _ in 1:batchsize];
outputs_gpu_3=[CUDA.zeros(n_select3,1) for _ in 1:batchsize];
svds_4=compute_single_svd_series(trial_matrices_trunc,inputs_gpu_3, outputs_gpu_3, 150,1, trials, batchsize)
svds_scaled_4=((svds_4) .^2 .-4 *(n))./(2(2(n)).^(1/3))

trial_matrices_trunc=[bidiag_trunc(n,n,n_select4) for _ in 1:trials]
svds_5=compute_single_svd_series(trial_matrices_trunc,inputs_gpu_3, outputs_gpu_3, 150,1, trials, batchsize)
svds_scaled_5=((svds_5) .^2 .-4 *(n))./(2(2(n)).^(1/3))


function tracy_wisdom!(dq,q,p, t) 
    dq[1]=q[2]
    dq[2]=t.*q[1]+(q[1].^3).*2
    dq[3]=q[4]
    dq[4]=q[1].^2
end
#source https://github.com/mitmath/18338/blob/master/JuliaNotebooks/TracyWidomTwoWays.ipynb
t0,tn = 5.0,-8.0
u0 = [airy(t0), airy(1, t0), 0, airy(t0)^2]
prob = ODEProblem(tracy_wisdom!, u0, (t0,tn))
sol = solve(prob,Tsit5(), reltol=1e-14, abstol=1e-14)
f(t,I,I′) = (t, -I′*exp(-I))

TracyWidomPDF_via_ODE(t) = f(t, sol(t)[[3,4]]...)[2]
plot(size=(800,400))
plot!( -5:.1:3, t -> TracyWidomPDF_via_ODE((t)), label="Tracy-Widom law")
plot!(svds_scaled_n, seriestype=:stephist, nbins=40,fill=true, fillalpha=0.5, normalize=true, 
label="p=100%")
plot!(svds_scaled_2, seriestype=:stephist, nbins=40,fill=true, fillalpha=0.5, normalize=true, 
label="p=5%")
plot!(svds_scaled_3, seriestype=:stephist, nbins=40,fill=true, fillalpha=0.5, normalize=true, 
label="p=2,5%")
plot!(svds_scaled_4, seriestype=:stephist, nbins=40,fill=true, fillalpha=0.5, normalize=true, 
label="p=1%", legend=:topleft)
plot!(svds_scaled_5, seriestype=:stephist, nbins=40,fill=true, fillalpha=0.5, normalize=true, 
label="p=0,5%", legend=:topleft)
plot!(title="Normalized largest singular values of bidiagonal matrix \n of size m=n="*string(n)* " for 150 Monte-Carlo trials")
savefig("limitTW.png")



################################################################################################
#-------------------------Random algorithms testing-----------------------------
################################################################################################

using MatrixMarket

M = Array(MatrixMarket.mmread("mbeacxc.mtx"))

M_vec=[M[:,i] for i in 1:size(M,2)]
M_norm=norm.(M_vec)
M_sortedi=sortperm(M_norm, rev=true)

k_values=unique(round.(Int,[10^i for i=0:0.2:2]))
M_S=svdvals((M))[k_values]
p_values=[10^i for i=0:-0.2:-2]
result_error=zeros(length(k_values),length(p_values))

for (i,p) in enumerate(p_values)
    println(p)
    n=round(Int,p*size(M,2))
    M_trunc=(M[:,M_sortedi[1:n]])
    M_truncvec=[M_trunc[i,:] for i in 1:size(M_trunc,1)]
    M_truncnorm=norm.(M_truncvec)
    M_sorteditrunc=sortperm(M_truncnorm, rev=true)
    M_trunc=M_trunc[M_sorteditrunc[1:n],:]

    M_trunc_S=svdvals(M_trunc)[k_values[k_values.<=n]]
    result_error[1:sum(k_values.<=n),i]=max.(abs.((M_S[k_values.<=n].-M_trunc_S)./M_S[k_values.<=n])',1e-7)
end

result_error=log10.(result_error)
cticks=[-2,-4,-6]
result_error=result_error[:,end:-1:1]
pyplot()
fontsize=10

plot(size=(600,400))
myplot=plot!(k_values, reverse(p_values), result_error', 
    xaxis=:log10, yaxis=:log10, 
    seriestype=:heatmap,  c=:viridis,
    fillalpha=0.5,     legend=:none, 
    xtickfontsize=fontsize,ytickfontsize=fontsize,
    ylabel="p",xlabel="k", yguidefontrotation=-90, 
    yguidefontsize=fontsize, xguidefontsize=fontsize,
    titlefontsize=fontsize,
    colorbar_tickfontsize=fontsize-1, 
    colorbar_ticks=(cticks, string.(round.(10. .^cticks,digits=6))),
    xlims=(1,100), ylims=(0.01,1)
    )

    savefig("examplersvd.png")