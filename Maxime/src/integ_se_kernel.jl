import GaussianProcesses: cov, cov!, get_params, EmptyData, Kernel, MatF64
import GaussianProcesses: set_params!, num_params, dKij_dθp

#=================
TYPE DEFINITIONS
=================#

type IntegSE <: Kernel
    ℓ2::Float64
    σ2::Float64
    IntegSE(ll::Float64, lσ::Float64) = new(exp(2*ll), exp(2*lσ))
end
IntegSE(se::SEIso) = IntegSE(get_params(se)...)

#========================
COVARIANCE AND DERIVATIVE
========================#

function cov(igse::IntegSE, x₁::Float64, x₂::Float64)
    u = x₁ / √(2 * igse.ℓ2)
    v = x₂ / √(2 * igse.ℓ2)
    uv = abs(u-v)
    ku = u * erf(u) + 1/√π * exp(-u^2)
    kv = v * erf(v) + 1/√π * exp(-v^2)
    kuv = uv * erf(uv) + 1/√π * exp(-uv^2)
    return igse.σ2 * igse.ℓ2 * √π * (ku + kv - kuv - 1/√π)
end
@inline function dk_dll(igse::IntegSE, x₁::Float64, x₂::Float64)
    u = x₁ / √(2 * igse.ℓ2)
    v = x₂ / √(2 * igse.ℓ2)
    uv = abs(u-v)
    dku = erf(u)
    dkv = erf(v)
    dkuv = erf(uv)
    
    ku = u * erf(u) + 1/√π * exp(-u^2)
    kv = v * erf(v) + 1/√π * exp(-v^2)
    kuv = uv * erf(uv) + 1/√π * exp(-uv^2)
    #dudll = -u
    dk_dll = igse.σ2 * igse.ℓ2 * √π * (-u*dku -v*dkv +uv*dkuv) +
             2 * cov(igse, x₁, x₂)
    return dk_dll
end

#========================
ALL OTHER METHODS
========================#
dk_dlσ(k::IntegSE, x₁::Float64, x₂::Float64) = 2.0*cov(k,x₁,x₂)

function cov(igse::IntegSE, x₁::AbstractVector{Float64}, x₂::AbstractVector{Float64})
    return cov(igse, x₁[1], x₂[1])
end
function cov!{M<:MatF64}(cK::MatF64, igse::IntegSE, X::M, data::EmptyData)
    dim,nobsv = size(X)
    @assert dim==1
    for i in 1:nobsv
        for j in 1:i
            cK[i,j] = cov(igse, X[1,i], X[1,j])
            cK[j,i] = cK[i,j]
        end
    end
end
function addcov!{M<:MatF64}(cK::MatF64, igse::IntegSE, X::M, data::EmptyData)
    dim,nobsv = size(X)
    @assert dim==1
    for i in 1:nobsv
        for j in 1:i
            cK[i,j] += cov(igse, X[1,i], X[1,j])
            cK[j,i] = cK[i,j]
        end
    end
end
function multcov!{M<:MatF64}(cK::MatF64, igse::IntegSE, X::M, data::EmptyData)
    dim,nobsv = size(X)
    @assert dim==1
    for i in 1:nobsv
        for j in 1:i
            cK[i,j] *= cov(igse, X[1,i], X[1,j])
            cK[j,i] = cK[i,j]
        end
    end
end
function cov{M<:MatF64}(igse::IntegSE, X::M, data::EmptyData)
    dim,nobsv = size(X)
    cK = Array(Float64, nobsv, nobsv)
    cov!(cK, igse, X, data)
    return cK
end

function cov!{M1<:MatF64,M2<:MatF64}(cK::MatF64, igse::IntegSE, X₁::M1, X₂::M2)
    @assert size(X₁,1)==size(X₂,1)==1
    nobsv1 = size(X₁, 2)
    nobsv2 = size(X₂, 2)
    for i in 1:nobsv1
        for j in 1:nobsv2
            cK[i,j] = cov(igse, X₁[1,i], X₂[1,j])
        end
    end
end  
function addcov!{M1<:MatF64,M2<:MatF64}(cK::MatF64, igse::IntegSE, X₁::M1, X₂::M2)
    @assert size(X₁,1)==size(X₂,1)==1
    nobsv1 = size(X₁, 2)
    nobsv2 = size(X₂, 2)
    for i in 1:nobsv1
        for j in 1:nobsv2
            cK[i,j] += cov(igse, X₁[1,i], X₂[1,j])
        end
    end
end  
function multcov!{M1<:MatF64,M2<:MatF64}(cK::MatF64, igse::IntegSE, X₁::M1, X₂::M2)
    @assert size(X₁,1)==size(X₂,1)==1
    nobsv1 = size(X₁, 2)
    nobsv2 = size(X₂, 2)
    for i in 1:nobsv1
        for j in 1:nobsv2
            cK[i,j] *= cov(igse, X₁[1,i], X₂[1,j])
        end
    end
end  
function cov{M1<:MatF64,M2<:MatF64}(igse::IntegSE, X₁::M1, X₂::M2)
    nobsv1 = size(X₁, 2)
    nobsv2 = size(X₂, 2)
    cK = Array(Float64, nobsv1, nobsv2)
    cov!(cK, igse, X₁, X₂)
    return cK
end

get_params(igse::IntegSE) = Float64[log(igse.ℓ2)/2.0, log(igse.σ2)/2.0]
get_param_names(::IntegSE) = [:ll, :lσ]
function set_params!(igse::IntegSE, hyp::Vector{Float64})
    length(hyp) == 2 || throw(ArgumentError("Squared exponential only has two parameters"))
    igse.ℓ2, igse.σ2 = exp(2.0*hyp)
end
num_params(igse::IntegSE) = 2

@inline function dKij_dθp{M<:MatF64}(igse::IntegSE, X::M, i::Int, j::Int, p::Int, dim::Int)
    if p==1
        return dk_dll(igse, X[1,i], X[1,j])
    elseif p==2
        return dk_dlσ(igse, X[1,i], X[1,j])
    else
        return NaN
    end
end
@inline function dKij_dθp{M<:MatF64}(igse::IntegSE, X::M, data::EmptyData, i::Int, j::Int, p::Int, dim::Int)
    return dKij_dθp(igse,X,i,j,p,dim)
end
kernel_data_key(::IntegSE, X::MatF64) = "EmptyData"
