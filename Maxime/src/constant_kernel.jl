using GaussianProcesses: Isotropic, MatF64, VecF64, SqEuclidean
import GaussianProcesses: cov, cov!, addcov!, kernel_data_key, get_params, set_params!, get_param_names, num_params, dk_dθp, dk_dlσ

type ConstantKernel <: Isotropic{SqEuclidean}
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
# cov{V1<:VecF64,V2<:VecF64}(k::ConstantKernel, x::V1, y::V2) = k.σ2

@inline function dk_dθp(k::ConstantKernel, r::Float64, p::Int)
    if p==1
        return dk_dlσ(se, r)
    else
        return NaN
    end
end
