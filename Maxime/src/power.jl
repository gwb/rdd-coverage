import Distributions: cdf, ccdf

function sim_power!(gpT::GP, gpC::GP, gpNull::GP, τ::Float64, 
            treat::BitVector, X∂::MatF64,
                chi_null::Vector{Float64}, 
                mLL_null::Vector{Float64}, 
                pval_invvar_null::Vector{Float64}
                ;
                update_mean=true)
    n = gpNull.nobsv
    null = MultivariateNormal(zeros(n), gpNull.cK)
    Ysim = rand(null)
    Ysim[treat] .+= τ

    gpT.y = Ysim[treat]
    gpC.y = Ysim[!treat]
    gpNull.y = Ysim
    if update_mean
        gpT.m = MeanConst(mean(gpT.y))
        gpC.m = MeanConst(mean(gpC.y))
        gpNull.m = MeanConst(mean(gpNull.y))
    end

    update_alpha!(gpT)
    update_alpha!(gpC)
    update_alpha!(gpNull)

    gpT.mLL = mLL(gpT)
    gpC.mLL = mLL(gpC)
    gpNull.mLL = mLL(gpNull)
    
    # mLL
    ΔmLL = gpT.mLL + gpC.mLL - gpNull.mLL
    pval_mLL = mean(mLL_null .> ΔmLL)
    
    # χ2
    pred_T = _predict_raw(gpT, gpT.k, gpT.m, X∂)
    pred_C = _predict_raw(gpC, gpC.k, gpC.m, X∂)
    μ = pred_T[1].-pred_C[1]
    Σ = pred_T[2]+pred_C[2]
    χ2 = dot(μ, Σ\μ)
    pval_χ2 = mean(chi_null .> χ2)
    
    # inverse-var
    pval_invvar_obs = pval_invvar(gpT, gpC, X∂)

    invvar_bootcalib = mean(pval_invvar_null .< pval_invvar_obs)
    invvar_calib = pval_invvar_calib(gpT, gpC, X∂)


    return (pval_mLL,pval_χ2,pval_invvar_obs,invvar_bootcalib,invvar_calib)
end

function nsim_power(gpT::GP, gpC::GP, τ::Float64, 
                X∂::MatF64,
                chi_null::Vector{Float64}, 
                mLL_null::Vector{Float64}, 
                pval_invvar_null::Vector{Float64},
                nsim::Int
                ;
                update_mean::Bool=true)
    gpT_mod = modifiable(gpT)
    gpC_mod = modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    power_sims = [sim_power!(gpT_mod, gpC_mod, gpNull, τ, treat, X∂, 
                            chi_null, mLL_null, pval_invvar_null
                            ;update_mean=update_mean) 
                  for _ in 1:nsim];
    return power_sims
end
