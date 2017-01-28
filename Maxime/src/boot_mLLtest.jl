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
immutable Point
    x::Float64
    y::Float64
end
function isLeft(a::Point, b::Point, c::Point)
     return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) > 0
end
function left_points(angle::Float64, shift::Float64, X::Matrix{Float64})
    dydx=tand(angle)
    perp=tand(angle+90)
    meanx = mean(X[1,:])
    meany = mean(X[2,:])
    direction=sign(cosd(angle+90))
    if direction==0.0
        direction=1.0
    end
    shift_x = shift*cosd(angle+90)*direction
    shift_y = shift*sind(angle+90)*direction
    a = Point(meanx+shift_x,meany+shift_y)
    if abs(dydx)<1.0
        dx=1.0
        dy=dydx
    else
        dx=1.0/dydx
        dy=1.0
    end
    b = Point(meanx+shift_x+dx,meany+dy+shift_y)
    
    n = size(X,2)
    points = [Point(X[1,i],X[2,i]) for i in 1:n]
    are_left = [isLeft(a,b,p) for p in points]
end

function plot_line(angle::Float64, shift::Float64, X::Matrix; kwargs...)
    meanx=mean(X[1,:])
    meany=mean(X[2,:])
    dydx=tand(angle)
    direction=sign(cosd(angle+90))
    if direction==0.0
        direction=1.0
    end
    shift_x = shift*cosd(angle+90)*direction
    shift_y = shift*sind(angle+90)*direction
    
    xlim=plt.xlim()
    ylim=plt.ylim()
    xlim_arr = Float64[xlim[1],xlim[2]]
    ylim_arr = Float64[ylim[1],ylim[2]]
    if dydx > 1e3
        plt.axvline(meanx+shift_x, color="red"; kwargs...)
    elseif dydx > 10
        plt.plot(meanx+(ylim_arr.-meany.-shift_y)./dydx+shift_x,ylim_arr; color="red", kwargs...)
    else
        plt.plot(xlim_arr,meany+(xlim_arr.-shift_x.-meanx).*dydx+shift_y; color="red", kwargs...)
    end
    plt.xlim(xlim)
    plt.ylim(ylim)
end

function prop_left(angle::Float64, shift::Float64, X::Matrix{Float64})
    are_left = left_points(angle, shift, X)
    return mean(are_left)
end

#copy-pasted from http://blog.mmast.net/bisection-method-julia
function bisection(f::Function, a::Number, b::Number;
                   tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a+b)/2
        fc = f(c)
        if fc == 0
            break
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return c
end

function shift_for_even_split(angle::Float64, X::Matrix)
    minx,maxx = extrema(X[1,:])
    return bisection(shift -> prop_left(angle,shift,X)-0.5, minx,maxx)
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

function placebo(angle::Float64, X::Matrix, Y::Vector, 
                 kern::Kernel, logNoise::Float64, 
                 nsim::Int; update_mean::Bool=false)
    shift = shift_for_even_split(angle, X)
    are_left = left_points(angle, shift, X)
    gp_left  = GP(X[:,are_left],  Y[are_left],  MeanConst(mean(Y[are_left])),  kern, logNoise)
    gp_right = GP(X[:,!are_left], Y[!are_left], MeanConst(mean(Y[!are_left])), kern, logNoise)
    pval = boot_mLLtest(gp_left, gp_right, nsim; update_mean=update_mean)
    return pval
end
