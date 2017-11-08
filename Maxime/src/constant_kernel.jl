using GaussianProcesses: Stationary, MatF64, VecF64, StationaryData, EmptyData
import GaussianProcesses: cov, cov!, addcov!, kernel_data_key, get_params, set_params!, get_param_names, num_params

type ConstantKernel <: Stationary
    σ2::Float64
    priors::Array
    ConstantKernel(lσ::Float64) = new(exp(2*lσ), [])
end
function set_params!(k::ConstantKernel, hyp::Vector{Float64})
    k.σ2 = exp(2*hyp[1])
end
get_params(k::ConstantKernel) = [0.5*log(k.σ2)]
get_param_names(::ConstantKernel) = [:lσ]
num_params(k::ConstantKernel) = 1
kernel_data_key{M<:MatF64}(k::ConstantKernel, X::M) = "ConstData"
cov(k::ConstantKernel, r::Float64) = k.σ2
cov{V1<:VecF64,V2<:VecF64}(k::ConstantKernel, x::V1, y::V2) = k.σ2

@inline function dk_dθp(k::ConstantKernel, r::Float64, p::Int)
    if p==1
        return dk_dlσ(se, r)
    else
        return NaN
    end
end
# distance{M<:MatF64}(k::ConstantKernel, X::M) = ones(Float64, size(X,2), size(X,2))
# distance{M1<:MatF64,M2<:MatF64}(k::ConstantKernel, X::M1, Y::M2) = ones(Float64, size(X,2), size(Y,2))
# distance{V1<:VecF64,V2<:VecF64}(k::ConstantKernel, x::V1, y::V2) = 1.0
function cov!{M<:MatF64}(cK::MatF64, k::ConstantKernel, X::M, data::EmptyData)
    fill!(cK, k.σ2)
end
function addcov!{M<:MatF64}(cK::MatF64, k::ConstantKernel, X::M, data::EmptyData)
    cK[:,:] += k.σ2
end   
