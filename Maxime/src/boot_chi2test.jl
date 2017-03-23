function _predict_raw{M<:MatF64,K<:Kernel,MF<:Mean}(gp::GP, k::K, m::MF, X::M)
    n = size(X, 2)
    cK = cov(k, X, gp.X)
    Lck = PDMats.whiten(gp.cK, cK')
    mu = mean(m,X) + cK*gp.alpha        # Predictive mean
    Sigma_raw = cov(k, X) - Lck'Lck     # Predictive covariance
    return mu, Sigma_raw
end
_predict_raw{M<:MatF64}(gp::GP, X::M) = _predict_raw(gp, gp.k, gp.m, X)

function chistat{K<:Kernel,MF<:Mean}(gpT::GP, gpC::GP, k::K, mT::MF, mC::MF, X∂::MatF64)
    pred_T = _predict_raw(gpT, k, mT, X∂)
    pred_C = _predict_raw(gpC, k, mC, X∂)
    μ = pred_T[1].-pred_C[1]
    Σ = pred_T[2]+pred_C[2]
    return dot(μ, Σ \ μ)
end
chistat(gpT::GP, gpC::GP, X∂::MatF64) = chistat(gpT,gpC,gpT.k,gpT.m,gpC.m,X∂)

function sim_chi_null{K<:Kernel,MF<:Mean}(gpT::GP, gpC::GP, gpNull::GP, k::K, mT::MF,mC::MF,
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

    return chistat(gpT, gpC, k, mT, mC, X∂)
end
function nsim_chi(gpT::GP, gpC::GP, X∂::MatF64, nsim::Int; update_mean::Bool=false)
    gpT_mod = modifiable(gpT)
    gpC_mod = modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GP([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)
    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    k = gpT_mod.k
    mT = gpT_mod.m
    mC = gpC_mod.m
    t_sims = [sim_chi_null(gpT_mod, gpC_mod, gpNull, k, mT,mC, treat, X∂; update_mean=update_mean) 
        for _ in 1:nsim];
    return t_sims
end

function boot_chi2test(gpT::GP, gpC::GP, X∂::MatF64, nsim::Int; update_mean::Bool=false)
    chi_obs = chistat(gpT, gpC, gpT.k, gpT.m, gpC.m, X∂)
    chi_sims = GeoRDD.nsim_chi(gpT, gpC, X∂, nsim; update_mean=update_mean)
    return mean(chi_sims .> chi_obs)
end


function placebo_chi(angle::Float64, X::MatF64, Y::Vector,
                 kern::Kernel, logNoise::Float64, 
                 nsim::Int; update_mean::Bool=false)
    shift = shift_for_even_split(angle, X)
    left = left_points(angle, shift, X)
    gp_left  = GP(X[:,left],  Y[left],  MeanConst(mean(Y[left])),  kern, logNoise)
    gp_right = GP(X[:,!left], Y[!left], MeanConst(mean(Y[!left])), kern, logNoise)
    X∂ = placebo_sentinels(angle, shift, X, 100)
    pval = boot_chi2test(gp_left, gp_right, X∂, nsim; update_mean=update_mean)
    return pval
end
