using RCall

function triangular_kernel(x1,x2,bw)
    max(1-abs(x1-x2)/bw, 0.0)
end

function IR_bandwidth(x::Vector, y::Vector, polydeg::Int)
    R"library(rdrobust)"
    bws_R = R"""rdbwselect($y,$x,p=$polydeg,bwselect="mserd")$bws"""
    bws = Vector(bws_R)
    bandwidth = bws[1]
    return bandwidth
end

function llr(x::Vector, y::Vector, polydeg::Int, xfit::Float64, bandwidth::Float64)
    n = length(y)
    X = Matrix{Float64}(n, polydeg+1)
    for deg in 0:polydeg
        X[:,deg+1] = x.^deg
    end
    w = triangular_kernel.(x, xfit, bandwidth)
    Xfit = [xfit^deg for deg in 0:polydeg]
    XtWX = X' * (w .* X)
    weights = (w .* X) * (XtWX \ Xfit)
    mu_fit = dot(weights, y)
    return weights, mu_fit
end

mutable struct LLR_RDD_Output
    x::Vector
    y::Vector
    polydeg::Int
    bandwidth::Float64
    τhat::Float64
    weights::Vector
end

function llr_rdd(x::Vector, y::Vector, polydeg::Int)
    n = length(y)
    bandwidth = IR_bandwidth(x, y, polydeg)
    treat = x .> 0
    ctrol = .!treat
    w_T, mu_T = llr(x[treat], y[treat], polydeg, 0.0, bandwidth)
    w_C, mu_C = llr(x[ctrol], y[ctrol], polydeg, 0.0, bandwidth)
    τhat = mu_T - mu_C
    weights = Vector{Float64}(n)
    weights[treat] = w_T
    weights[ctrol] = -w_C
    return LLR_RDD_Output(
        x, y, polydeg,
        bandwidth,
        τhat,
        weights
        )
end
