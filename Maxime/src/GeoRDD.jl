module GeoRDD
    using GaussianProcesses
    using PDMats
    using Optim
    include("multigp_covars.jl")
    include("GPrealisations.jl")
    include("point_estimates.jl")
    include("cliff_face.jl")
end
