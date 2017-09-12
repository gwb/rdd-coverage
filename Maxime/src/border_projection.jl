import GeoInterface: coordinates
import LibGEOS
import IterTools
import DataStructures
import LibGEOS: nearestPoints, interpolate, distance
import GaussianProcesses: GPE

function projection_points(gp::GPE, border::BorderType; maxdist::Float64=Inf)
    X∂ = Array{Float64}(2, gp.nobsv)
    distances = Vector{Float64}(gp.nobsv)
    for i in 1:gp.nobsv
        # obtain coordinates for treatment point
        x,y = gp.X[:,i]
        point = LibGEOS.Point(x,y)
        # projection onto border (as distance along border)
        distances[i] = distance(border, point)
        proj_point = nearestPoints(border, point)[1]
        # get border point from distance
        # proj_point = interpolate(border, proj_dist)
        # get coordinates from point
        proj_x, proj_y = coordinates(proj_point)
        X∂[1,i] = proj_x
        X∂[2,i] = proj_y
    end
    return X∂[:, distances .<= maxdist]
end

function proj_estimator(gpT::GPE, gpC::GPE, border::BorderType)
    X∂_treat = projection_points(gpT, border)
    X∂_ctrol = projection_points(gpC, border)
    X∂ = [X∂_treat X∂_ctrol]
    
    μpost, Σpost = cliff_face(gpT, gpC, X∂)
    return unweighted_mean(μpost, Σpost)
end

function buffer_infinite_proj_sentinels(gpT::GPE, gpC::GPE, border::BorderType,
                                    buffer::Float64, gridspace::Float64;
                                    density = ((s1,s2) -> 1.0))
    # Add a buffer around the border, returns a "border area".
    #=buffer_polygon = LibGEOS.buffer(border, buffer)=#

    # Obtain a convex hull containing all the data.
    X_multi_treat = LibGEOS.MultiPoint([gpT.X[:, i] for i in 1:gpT.nobsv])
    X_multi_ctrol = LibGEOS.MultiPoint([gpC.X[:, i] for i in 1:gpC.nobsv])
    convexhull_treat = LibGEOS.convexhull(X_multi_treat)
    convexhull_ctrol = LibGEOS.convexhull(X_multi_ctrol)
    data_hull = LibGEOS.union(convexhull_treat, convexhull_ctrol)

    # Obtain grid of points in bounding box for all data.
    X1_min = min(minimum(gpT.X[1,:]), minimum(gpC.X[1,:]))
    X1_max = max(maximum(gpT.X[1,:]), maximum(gpC.X[1,:]))
    X1_grid = X1_min:gridspace:X1_max

    X2_min = min(minimum(gpT.X[2,:]), minimum(gpC.X[2,:]))
    X2_max = max(maximum(gpT.X[2,:]), maximum(gpC.X[2,:]))
    X2_grid = X2_min:gridspace:X2_max

    projected_weights = Dict{Vector{Float64}, Float64}()
    for (s1,s2) in IterTools.product(X1_grid, X2_grid)
        # convert to a point obejct
        p = LibGEOS.Point(s1,s2)
        # Only keep points that are both within `buffer` of
        # the border, and are in the convex hull of the data.
        if LibGEOS.distance(p, border) > buffer
            continue
        elseif !LibGEOS.within(p, data_hull)
            continue
        end
        # compute the population density at this location
        ρ = density(s1,s2)

        # project the point onto the border
        projected = LibGEOS.nearestPoints(p, border)[2]
        # obtain the coordinates of the projected point
        proj_coords = GeoInterface.coordinates(projected)[1:2]
        if haskey(projected_weights, proj_coords)
            projected_weights[proj_coords] += ρ
        else
            # initialize
            projected_weights[proj_coords] = ρ
        end
    end

    unique_projected_points = collect(keys(projected_weights))
    X∂_projected = [[p[1] for p in unique_projected_points]';
                    [p[2] for p in unique_projected_points]']
    # And the counts as the weights.
    weights = collect(values(projected_weights))
    # Note: assumes keys() and values() are in the same order
    return X∂_projected, weights
end

function buffer_infinite_proj_estim(gpT::GPE, gpC::GPE, border::BorderType,
                                    buffer::Float64, gridspace::Float64;
                                    density = ((s1,s2) -> 1.0))
    X∂_projected, weights = buffer_infinite_proj_sentinels(
            gpT, gpC, border, buffer, gridspace; density=density
            )
            
    # Obtain the cliff face for these sentinels...
    μpost, Σpost = cliff_face(gpT, gpC, X∂_projected)
    # ... and apply a weighted mean estimator.
    return weighted_mean(μpost, Σpost, weights)
end
