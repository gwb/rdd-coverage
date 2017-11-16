import Distributions: cdf, ccdf

function sim_power!(gpT::GPE, gpC::GPE, gpNull::GPE, τ::Float64, 
            treat::BitVector, Xb::MatF64,
                chi_null::Vector{Float64}, 
                mll_null::Vector{Float64}, 
                pval_invvar_null::Vector{Float64},
                Σcliff::PDMat,
                cK_T::MatF64,
                cK_C::MatF64,
                KCT::MatF64
                ;
                update_mean=true)
    # simulate data
    n = gpNull.nobsv
    null = MultivariateNormal(zeros(n), gpNull.cK)
    Ysim = rand(null)
    Ysim[treat] .+= τ

    # Modify data in GP objects
    gpT.y[:] = Ysim[treat]
    gpC.y[:] = Ysim[.!treat]
    gpNull.y[:] = Ysim
    if update_mean
        gpT.m = MeanConst(mean(gpT.y))
        gpC.m = MeanConst(mean(gpC.y))
        gpNull.m = MeanConst(mean(gpNull.y))
    end

    # Update GPs
    update_alpha!(gpT)
    update_alpha!(gpC)
    update_alpha!(gpNull)

    gpT.mll = mll(gpT)
    gpC.mll = mll(gpC)
    gpNull.mll = mll(gpNull)
    
    # mll
    Δmll = gpT.mll + gpC.mll - gpNull.mll
    pval_mll = mean(mll_null .> Δmll)
    
    # χ2
    χ2 = chistat(gpT, gpC, Xb, Σcliff, cK_T, cK_C)
    pval_χ2 = mean(chi_null .> χ2)
    
    # inverse-var
    pval_invvar_obs = pval_invvar(gpT, gpC, Xb, Σcliff, cK_T, cK_C)

    invvar_bootcalib = mean(pval_invvar_null .< pval_invvar_obs)
    invvar_calib = pval_invvar_calib(gpT, gpC, Xb, Σcliff, cK_T, cK_C, KCT)

    return (pval_mll,pval_χ2,pval_invvar_obs,invvar_bootcalib,invvar_calib)
end

function nsim_power(gpT::GPE, gpC::GPE, τ::Float64, 
                Xb::MatF64,
                chi_null::Vector{Float64}, 
                mll_null::Vector{Float64}, 
                pval_invvar_null::Vector{Float64},
                nsim::Int
                ;
                update_mean::Bool=true)
    gpT_mod = modifiable(gpT)
    gpC_mod = modifiable(gpC)
    yNull = [gpT.y; gpC.y]
    gpNull = GPE([gpT.X gpC.X], yNull, MeanConst(mean(yNull)), gpT.k, gpT.logNoise)

    _, Σcliff = cliff_face(gpT, gpC, Xb)
    cK_T = cov(gpT.k, Xb, gpT.X)
    cK_C = cov(gpC.k, Xb, gpC.X)
    KCT = cov(gpC.k, gpC.X, gpT.X)

    treat = BitVector(gpNull.nobsv)
    treat[:] = false
    treat[1:gpT.nobsv] = true
    power_sims = [sim_power!(gpT_mod, gpC_mod, gpNull, τ, treat, Xb, 
                            chi_null, mll_null, pval_invvar_null,
                            Σcliff, cK_T, cK_C, KCT,
                            ;update_mean=update_mean) 
                  for _ in 1:nsim];
    return power_sims
end
