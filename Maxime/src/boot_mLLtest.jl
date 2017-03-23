function mLL(gp::GP)
    μ = mean(gp.m,gp.X)
    return -dot((gp.y - μ),gp.alpha)/2.0 - logdet(gp.cK)/2.0 - gp.nobsv*log(2π)/2.0 # Marginal log-likelihood
end

function sim_logP(gpT::GP, gpC::GP, gpNull::GP, treat::BitVector; update_mean::Bool=false)
    n = gpNull.nobsv
    null = MultivariateNormal(zeros(n), gpNull.cK)
    Ysim = rand(null)
    
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
    
    mLL_altv = gpT.mLL + gpC.mLL
    mLL_null = gpNull.mLL
    return mLL_null, mLL_altv
end

function nsim_logP(gpT::GP, gpC::GP, 
                  nsim::Int; update_mean::Bool=false)
    gpT_mod = modifiable(gpT)
    gpC_mod = modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    mLL_sims = [sim_logP(gpT_mod, gpC_mod, gpNull, treat; update_mean=update_mean) 
        for _ in 1:nsim];
    return mLL_sims
end

function data_logP(gpT::GP, gpC::GP)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    update_mll!(gpNull)
    return gpNull.mLL, gpT.mLL+gpC.mLL
end

function boot_mLLtest(gpT::GP, gpC::GP, nsim::Int; update_mean::Bool=false)
    mLL_alt = gpT.mLL + gpC.mLL
    
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    
    mLL_null = gpNull.mLL
    mLL_sims = GeoRDD.nsim_logP(gpT, gpC, nsim; update_mean=update_mean)
    
    mLL_sim_null = [sim[1] for sim in mLL_sims]
    mLL_sim_altv = [sim[2] for sim in mLL_sims]
    
    ΔmLL_obs = mLL_alt - mLL_null
    ΔmLL_sim = mLL_sim_altv .- mLL_sim_null
    return mean(ΔmLL_sim .> ΔmLL_obs)
end

function placebo_mLL(angle::Float64, X::Matrix, Y::Vector, 
                 kern::Kernel, logNoise::Float64, 
                 nsim::Int; update_mean::Bool=false)
    shift = shift_for_even_split(angle, X)
    left = left_points(angle, shift, X)
    gp_left  = GP(X[:,left],  Y[left],  MeanConst(mean(Y[left])),  kern, logNoise)
    gp_right = GP(X[:,!left], Y[!left], MeanConst(mean(Y[!left])), kern, logNoise)
    pval = boot_mLLtest(gp_left, gp_right, nsim; update_mean=update_mean)
    return pval
end
