using Optim
using GaussianProcesses: Mean, Kernel, evaluate, metric
import GaussianProcesses: optimize!, get_optim_target
import GaussianProcesses: num_params, set_params!, get_params
import GaussianProcesses: update_mll_and_dmll!, update_mll!

type GPRealisations
    reals::Vector{GP}
    k::Kernel
    logNoise::Float64
    mLL::Float64
    dmLL::Vector{Float64}
end

function GPRealisations(reals::Vector{GP})
    first = reals[1]
    gpr = GPRealisations(reals, first.k, first.logNoise, NaN, [])
end

function get_params(gpr::GPRealisations; noise::Bool=true, mean::Bool=true, kern::Bool=true)
    params = Float64[]
    if noise; push!(params, gpr.logNoise); end
    if mean
        for gp in gpr.reals
            append!(params, get_params(gp.m))
        end
    end
    if kern; append!(params, get_params(gpr.k)); end
    return params
end
function propagate_params!(gpr::GPRealisations; noise::Bool=true, kern::Bool=true)
    for gp in gpr.reals
        # harmonize parameters
        if kern
            gp.k = gpr.k
        end
        if noise
            gp.logNoise = gpr.logNoise
        end
    end
end
function set_params!(gpr::GPRealisations, hyp::Vector{Float64}; noise::Bool=true, mean::Bool=true, kern::Bool=true)
    # println("mean=$(mean)")
    istart=1
    if noise
        gpr.logNoise = hyp[istart]
        istart += 1
    end
    if mean
        for gp in gpr.reals
            set_params!(gp.m, hyp[istart:istart+num_params(gp.m)-1])
            istart += num_params(gp.m)
        end
    end
    if kern
        set_params!(gpr.k, hyp[istart:end])
    end
    propagate_params!(gpr, noise=noise, kern=kern)
end

function update_mll!(gpr::GPRealisations)
    mLL = 0.0
    for gp in gpr.reals
        update_mll!(gp)
        mLL += gp.mLL
    end
    gpr.mLL = mLL
    return mLL
end
function update_mll_and_dmll!(gpr::GPRealisations, Kgrads::Dict{Int,Matrix}, ααinvcKIs::Dict{Int,Matrix}; 
                              noise::Bool=true, mean::Bool=true, kern::Bool=true)
    gpr.mLL = 0.0
    gpr.dmLL = zeros(get_params(gpr; noise=noise, mean=mean, kern=kern))
    imean=2
    ikern=collect((length(gpr.dmLL)-num_params(gpr.k)+1):length(gpr.dmLL))
    for gp in gpr.reals
        update_mll_and_dmll!(gp, Kgrads[gp.nobsv], ααinvcKIs[gp.nobsv]; 
            noise=noise,mean=mean,kern=kern)
        gpr.mLL += gp.mLL
        dmLL_indices=Int[]
        if noise
            push!(dmLL_indices, 1)
        end
        if mean
            append!(dmLL_indices, imean:imean+num_params(gp.m)-1)
            imean+=num_params(gp.m)
        end
        if kern
            append!(dmLL_indices, ikern)
        end
        gpr.dmLL[dmLL_indices] .+= gp.dmLL
    end
    return gpr.dmLL
end
function update_mll_and_dmll!(gpr::GPRealisations; kwargs...)
    Kgrads = Dict{Int,Matrix}()
    ααinvcKIs = Dict{Int,Matrix}()
    for gp in gpr.reals
        if haskey(Kgrads, gp.nobsv)
            continue
        end
        Kgrads[gp.nobsv] = Array(Float64, gp.nobsv, gp.nobsv)
        ααinvcKIs[gp.nobsv] = Array(Float64, gp.nobsv, gp.nobsv)
    end
    return update_mll_and_dmll!(gpr, Kgrads, ααinvcKIs)
end

function get_optim_target(gpr::GPRealisations; noise::Bool=true, mean::Bool=true, kern::Bool=true)
    Kgrads = Dict{Int,Matrix}()
    ααinvcKIs = Dict{Int,Matrix}()
    for gp in gpr.reals
        if haskey(Kgrads, gp.nobsv)
            continue
        end
        Kgrads[gp.nobsv] = Array(Float64, gp.nobsv, gp.nobsv)
        ααinvcKIs[gp.nobsv] = Array(Float64, gp.nobsv, gp.nobsv)
    end
    function mll(hyp::Vector{Float64})
        try
            set_params!(gpr, hyp; noise=noise, mean=mean, kern=kern)
            update_mll!(gpr)
            return -gpr.mLL
        catch err
             if !all(isfinite(hyp))
                println(err)
                return Inf
            elseif isa(err, ArgumentError)
                println(err)
                return Inf
            elseif isa(err, Base.LinAlg.PosDefException)
                println(err)
                return Inf
            else
                throw(err)
            end
        end        
    end

    function mll_and_dmll!(hyp::Vector{Float64}, grad::Vector{Float64})
        try
            set_params!(gpr, hyp; noise=noise, mean=mean, kern=kern)
            update_mll_and_dmll!(gpr, Kgrads, ααinvcKIs; noise=noise, mean=mean, kern=kern)
            grad[:] = -gpr.dmLL
            return -gpr.mLL
        catch err
             if !all(isfinite(hyp))
                println(err)
                return Inf
            elseif isa(err, ArgumentError)
                println(err)
                return Inf
            elseif isa(err, Base.LinAlg.PosDefException)
                println(err)
                return Inf
            else
                throw(err)
            end
        end 
    end
    function dmll!(hyp::Vector{Float64}, grad::Vector{Float64})
        mll_and_dmll!(hyp::Vector{Float64}, grad::Vector{Float64})
    end

    func = DifferentiableFunction(mll, dmll!, mll_and_dmll!)
    return func
end
function optimize!(gpr::GPRealisations; noise::Bool=true, mean::Bool=true, kern::Bool=true,
                    method=ConjugateGradient(), kwargs...)
    func = get_optim_target(gpr, noise=noise, mean=mean, kern=kern)
    init = get_params(gpr;  noise=noise, mean=mean, kern=kern)  # Initial hyperparameter values
    results=optimize(func,init; method=method, kwargs...)                     # Run optimizer
    set_params!(gpr, Optim.minimizer(results), noise=noise,mean=mean,kern=kern)
    update_mll!(gpr)
    return results
end
