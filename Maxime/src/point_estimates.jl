using Distributions: Normal
using PDMats: AbstractPDMat

function inverse_variance(μ::AbstractVector, Σ::AbstractMatrix)
    n = size(μ)
    denom = sum(Σ \ ones(n))
    τhat = sum(Σ \ μ) / denom
    Vτhat = 1.0/denom
    τpost=Normal(τhat, √Vτhat)
    return τpost
end
    
# It really should be possible to write this function only once,
# but for some reason AbstractPDMat does not derive from
# AbstractMatrix
function inverse_variance(μ::AbstractVector, Σ::AbstractPDMat)
    n = size(μ)
    denom = sum(Σ \ ones(n))
    τhat = sum(Σ \ μ) / denom
    Vτhat = 1.0/denom
    τpost=Normal(τhat, √Vτhat)
    return τpost
end
