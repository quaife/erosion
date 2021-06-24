using Parameters
# The set of parameters.
@with_kw struct ParamSet
	# The data file for the initial geometry
	infolder::AbstractString = "input_geos/test/"
	label::AbstractString = "02-1"
	# Varied computational parameters
	npts::Int = 128			# The number of points per body, default 128.
	ibary::Int = 1			# Use barycentric (1) or not (0).
	ifmm::Int = 1			# Use the FFM (1) or not (0).
	ibc::Int = 1			# Use slip BCs (1) or no-slip (0)
	# Varied physical parameters
	epsfac::Float64 = 15	# Smoothing parameter for curvature driven flow.
	sigfac::Float64 = 10	# Smoothing parameter for the stress.
	dt::Float64 = 1e-3		# The time step.
	outstride::Int = 4		# The number of steps per output, default 4.
	# Fixed physical parameters
	fixpdrop::Bool = 1		# Fix the pressure drop (1) or not (0)
	fixarea::Bool = 0		# Keep the area fixed (1) or not (0)
	tfin::Float64 = 1.0		# The final time.
	# Fixed computational parameters
	maxl::Int = 8000		# The maximum number of GMRES iterations.
	nouter::Int = 1024		# The number of points on the outer boundary, default 1024.
	# Derived parameters
	epsilon::Float64 = epsfac/npts
	sigma::Float64 = sigfac/npts
end

# To call the main routine:
#include("main.jl")
#main(ParamSet())