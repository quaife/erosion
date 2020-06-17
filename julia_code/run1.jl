using Parameters
include("main.jl")

@with_kw struct paramset
	# The data file for the initial geometry
	infolder::AbstractString = "input_geos/afrac06/"
	label::AbstractString = "20-5"
	# Varied computational parameters
	npts::Int = 128			# The number of points per body.
	ibary::Int = 1			# Use barycentric (1) or not (0).
	iffm::Int = 1			# Use the FFM (1) or not (0).
	ibc::Int = 1			# Use slip BCs (1) or no-slip (0)
	# Varied physical parameters
	epsfac::Float64 = 15	# Smoothing parameter for curvature driven flow.
	sigfac::Float64 = 10	# Smoothing parameter for the stress.
	dt::Float64 = 1e-3		# The time step.
	nout::Int = 4			# The number of steps per output.
	# Fixed physical parameters
	fixpdrop::Int = 1		# Fix the pressure drop (1) or not (0)
	fixarea::Int = 0		# Keep the area fixed (1) or not (0)
	tfin::Float64 = 1.0		# The final time.
	# Fixed computational parameters
	maxl::Int = 8000		# The maximum number of GMRES iterations.
	nouter::Int = 1024		# The number of points on the outer boundary.
	# Derived parameters
	infile::AbstractString = string("../",infolder,label,".circ")
	epsilon::Float64 = epsfac/npts
	sigma::Float64 = sigfac/npts
end

cputime = @elapsed( erosion(paramset()) )

cputime = sig(cputime/3600,3)
println("\n\n\nCOMPLETED SIMULATION")
println("cpu time = ", cputime, " hours.\n\n")

## Other things to keep track of in addition to fixed parameters: cputime, cntout, 


