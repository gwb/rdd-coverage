module GeoRDD
    using GaussianProcesses
    using PDMats
    using Optim
    include("constant_kernel.jl")
    include("multigp_covars.jl")
    include("GPrealisations.jl")
    include("point_estimates.jl")
    include("cliff_face.jl")
    include("placebo_geometry.jl")
    include("boot_chi2test.jl")
    include("boot_mLLtest.jl")
    include("boot_invvar_test.jl")
    include("power.jl")
end
