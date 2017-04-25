function pval_invvar(gpT::GP, gpC::GP, X∂::MatF64)
    μ, Σ = cliff_face(gpT, gpC, X∂)
    invvar = inverse_variance(μ, Σ)
    pval_invvar = 2*min(cdf(invvar, 0.0), ccdf(invvar, 0.0))
    return pval_invvar
end

function sim_invvar(gpT::GP, gpC::GP, gpNull::GP, 
            treat::BitVector, X∂::MatF64; update_mean::Bool=false)
    n = gpNull.nobsv
    null = MultivariateNormal(zeros(n), gpNull.cK)
    Ysim = rand(null)

    gpT.y = Ysim[treat]
    gpC.y = Ysim[!treat]

    if update_mean
        gpT.m = mT = MeanConst(mean(gpT.y))
        gpC.m = mC = MeanConst(mean(gpC.y))
    end

    update_alpha!(gpT)
    update_alpha!(gpC)

    return pval_invvar(gpT, gpC, X∂)
end

function nsim_invvar_pval(gpT::GP, gpC::GP, X∂::MatF64, nsim::Int; update_mean::Bool=false)
    gpT_mod = modifiable(gpT)
    gpC_mod = modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    pval_sims = Float64[sim_invvar(gpT_mod, gpC_mod, gpNull, treat, X∂;
        update_mean=update_mean) 
        for _ in 1:nsim];
    return pval_sims
end

function boot_invvar(gpT::GP, gpC::GP, X∂::MatF64, nsim::Int; update_mean::Bool=false)
    pval_obs = pval_invvar(gpT, gpC, X∂)
    pval_sims = sim_invvar(gpT, gpC, X∂, nsim; update_mean=update_mean)
    return mean(pval_obs .< pval_sims)
end

#============================================
    ANALYTIC INSTEAD OF BOOTSTRAP CALIBRATION
=============================================#
function pval_invvar_calib(gpT::GP, gpC::GP, X∂::Matrix)
    extrap◫_T = GaussianProcesses.predict(gpT, X∂; full_cov=true)
    extrap◫_C = GaussianProcesses.predict(gpC, X∂; full_cov=true)
    μ∂ = extrap◫_T[1].-extrap◫_C[1]
    n = size(μ∂)
    Σ∂ = extrap◫_T[2]+extrap◫_C[2]
    
    K∂C = cov(gpC.k, X∂, gpC.X)
    KC∂ = K∂C'
    KCC = gpC.cK
    
    KCT = cov(gpC.k, gpC.X, gpT.X)
    KTT = gpT.cK
    KT∂ = cov(gpT.k, gpT.X, X∂)
    K∂T = KT∂'

    AT_c = KTT \ KT∂
    AC_c = KCC \ KC∂
    AT = AT_c'
    AC = AC_c'
    cov_μδ = AT*full(KTT)*AT_c + AC*full(KCC)*AC_c - AC*full(KCT)*AT_c - AT*full(KCT)'*AC_c
    
    cov_μτ= sum((Σ∂ \ cov_μδ) * (Σ∂ \ ones(n)))
    null = Normal(0.0, √cov_μτ)
    
    μτ_numer = sum(Σ∂ \ μ∂) # numerator only
    
    pval = 2*ccdf(null, abs(μτ_numer))
    return pval
end
function sim_invvar_calib(gpT::GP, gpC::GP, gpNull::GP, 
            treat::BitVector, X∂::MatF64; update_mean::Bool=false)
    n = gpNull.nobsv
    null = MultivariateNormal(zeros(n), gpNull.cK)
    Ysim = rand(null)

    gpT.y = Ysim[treat]
    gpC.y = Ysim[!treat]

    if update_mean
        gpT.m = mT = MeanConst(mean(gpT.y))
        gpC.m = mC = MeanConst(mean(gpC.y))
    end

    GeoRDD.update_alpha!(gpT)
    GeoRDD.update_alpha!(gpC)

    return pval_invvar_calib(gpT, gpC, X∂)
end
function nsim_invvar_calib(gpT::GP, gpC::GP, X∂::MatF64, nsim::Int; update_mean::Bool=false)
    gpT_mod = GeoRDD.modifiable(gpT)
    gpC_mod = GeoRDD.modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    pval_sims = Float64[sim_invvar_calib(gpT_mod, gpC_mod, gpNull, treat, X∂;
        update_mean=update_mean) 
        for _ in 1:nsim];
    return pval_sims
end

function placebo_invvar(angle::Float64, X::MatF64, Y::Vector,
                 kern::Kernel, logNoise::Float64)
    shift = shift_for_even_split(angle, X)
    left = left_points(angle, shift, X)
    gp_left  = GP(X[:,left],  Y[left],  MeanConst(mean(Y[left])),  kern, logNoise)
    gp_right = GP(X[:,!left], Y[!left], MeanConst(mean(Y[!left])), kern, logNoise)
    X∂ = placebo_sentinels(angle, shift, X, 100)
    pval = pval_invvar_calib(gp_left, gp_right, X∂)
    return pval
end
