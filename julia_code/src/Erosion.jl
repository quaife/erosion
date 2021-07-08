# top-level of the package, ties Julia files and modules together.
module Erosion

    __precompile__(true)

    # load submodules
    include("SpectralMethods.jl") # FFTs
    using .SpectralMethods
    include("ThetaLen.jl") # data structures
    using .ThetaLen
    include("DensityStress.jl") # FORTRAN wrappers
    using .DensityStress
    include("TimeStepping.jl") # timestepping
    using .TimeStepping

    # global includes
    include("Main.jl") # entry point and global functions

    # load extension modules
    include("MakeGeos.jl") # generate geometry
    using .MakeGeos
    include("Postprocessing.jl")
    using .Postprocessing

    # expose functionality to user
    export run_erosion, ParamSet, ThLenDenType, ThetaLenType

end
