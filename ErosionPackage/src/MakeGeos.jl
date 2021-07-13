#= OBJECTIVE: Generate initial geometries. =#
# To run: makegeos(10, 0.5, 4)
# To make multiple geometries run: make9geos(10, 0.5)

module MakeGeos
export make_geos, make9geos, CircType, plotcircs, get_afrac

using Erosion.ThetaLen: getalpha
using Random
using Distributions
using LinearAlgebra: norm
using Plots
using JLD2

#--------------- INITIALIZATION ---------------#
geosfolder() = "input_geos/"
figfolder() = "zFigsGeos/"

mutable struct CircType
	rad::Float64; xc::Float64; yc::Float64
end
#-------------------------------------------------#

#--------------- SUPPORTING ROUTINES ---------------#
# The repulsive force on circ1 due to circ2.
function fcircs(circ1::CircType, circ2::CircType, pow::Float64, buff::Float64, )
	# The critical distance.
	dcrit = (1+buff)*(circ1.rad + circ2.rad)
	# The vector from circ2 center to circ1 center, and the distance.
	v21 = [circ1.xc - circ2.xc, circ1.yc - circ2.yc]
	dist = norm(v21)
	# Compute the magnitude of the repulsive force if the bodies are within distance dcrit.
	dist < dcrit ? fmag = (1-dist/dcrit)^pow : fmag = 0.
	# Compute the repulsive force; if the bodies are right on top of each other, set the force to zero.
	f21 = fmag*v21/dist
	dist < 100*eps(dcrit) ? f21 *= 0 : 0
	return f21
end

# The repulsive force due to the four walls.
function fwall(circ::CircType, pow::Float64, buff::Float64)
	rcrit = (1+4*buff)*circ.rad
	# Get the force component if the x or y coordinate is within rcrit if any wall.
	function fcomp(xx::Float64)
		ff = 0.
		if xx < -1+rcrit
			ff = 1 - (xx+1)/rcrit
		elseif xx > 1-rcrit
			ff = -1 + (1-xx)/rcrit
		end
		return abs(ff)^pow * sign(ff)
	end
	# Return the force component in each direction.
	return [fcomp(circ.xc), fcomp(circ.yc)]
end

# The sum of all forces on each circle.
function forcesum(circvec::Vector{CircType}, pow::Float64, buff::Float64)
	nbods = length(circvec)
	fx, fy = [zeros(Float64, nbods) for nn=1:2]
	# Compute the total forces.
	for nn = 1:nbods
		fn = [0., 0.]
		# Compute the forces due to the other circles.
		for mm = 1:nbods
			mm != nn ? fn += fcircs(circvec[nn], circvec[mm], pow, buff) : 0
		end
		# Compute the forces due to the walls.
		fn += fwall(circvec[nn], pow, buff)
		fx[nn], fy[nn] = fn[1], fn[2]
	end
	return fx, fy
end

# Shift the circles with given force.
function shiftcircs(circvec::Vector{CircType}, fxv::Vector{Float64}, fyv::Vector{Float64}, dt::Float64, sigma::Float64)
	nbods = length(circvec)
	rvec = randn(2, nbods)
	rfac = sigma*sqrt(dt)
	for bod = 1:nbods
		circvec[bod].xc += dt*fxv[bod] + rfac*rvec[1, bod]
		circvec[bod].yc += dt*fyv[bod] + rfac*rvec[2, bod]
	end
end

# Plot the circles.
function plotcircs(circvec::Vector{CircType}, nfile::Integer, seed::Integer)
	npts = 128
	width = 400; height = 400
	nbods = length(circvec)
	alpha = getalpha(npts)
	# Set name of the folder and file.
	if nfile >= 0
		figname = string(figfolder(), "circ", lpad(nfile,4,"0"), ".pdf")
	else
		figname = string(geosfolder(), lpad(nbods,2,"0"), "-", seed, ".pdf")
	end
	# Make the figure.
	pp = plot(xlim=(-1,1), ylim=(-1,1), size=(width,height), leg=false);
	for bod = 1:nbods
		circ = circvec[bod]
		xx = circ.xc .+ circ.rad*cos.(alpha)
		yy = circ.yc .+ circ.rad*sin.(alpha)
		plot!(pp, xx, yy, color="black");
	end
	savefig(pp, figname)
end

# Compute the area fraction of the circles.
get_afrac(radvec::Vector{<:AbstractFloat}) = 0.25*pi*sum(radvec.^2)
# Dispatch for vector of CircType.
function get_afrac(circvec::Vector{CircType})
	radvec = [circvec[nn].rad for nn = 1:length(circvec)]
	return get_afrac(radvec)
end
#-------------------------------------------------#

#--------------- MAIN ROUTINES ---------------#
# Main routine to make the geometry.
function make_geos(nbods::Int, areafrac::Float64, seed::Int=1; makeplots::Bool = true)
	#= Parameters: 
	buff sets the buffer as a percentage of distance between bodies; beyond that distance there is no force.
	bolap sets an analogous buffer to test if the bodies are overlapping.
	pow sets the power-law behavior of the repulsive force between bodies.
	dt sets the time step in the simulated annealing. =#
	buff = 0.08
	bolap = 0.03
	pow = 0.5
	dt = 5e-2

	# Check that the desired area fraction is not too high; 0.91 is the absolute upper bound. 
	@assert areafrac < 0.71
	# Seed the random number generator and create the list of random radii.
	println("\n\nseed = ", seed)	
	Random.seed!(seed)
	# Two options for the radius distribution are Rayleigh and Chi.
	dray = Rayleigh(); dchi = Chi(4); drad = dchi
	radvec = rand(drad, nbods)
	# Rescale the radii to achieve desired area fraction.
	radvec *= sqrt( areafrac / get_afrac(radvec) )
	# Chose the provisional centers from a uniform distribution.
	duni = Uniform(-1,1)
	xc, yc = rand(duni, nbods), rand(duni, nbods)

	# Create the list of circles
	circvec = [CircType(radvec[bod], xc[bod], yc[bod]) for bod=1:nbods]
	# Create the folder for plots.
	if isdir(figfolder()) rm(figfolder(); recursive=true) end; mkdir(figfolder())

	# Shift the centers until no overlap.
	cnt = 0; pass = true; foverlap = 1.0
	while(cnt < 50 || foverlap > 1e-10)
		# Plot the circles.
		if mod(cnt, 10) == 0
			println("count = ", cnt)
			if makeplots plotcircs(circvec, cnt, seed) end
		end
		# Compute the forces and shift the circles.
		fx, fy = forcesum(circvec, pow, buff)
		sigma = 0.15*exp(-0.5*cnt*dt)
		shiftcircs(circvec, fx, fy, dt, sigma)
		# Test if the bodies are nearly overlapping with the overlap buffer.
		fxo, fyo = forcesum(circvec, pow, bolap)
		foverlap = max(norm(fxo, Inf), norm(fyo, Inf))
		# println("foverlap = ", round(foverlap, sigdigits=3))
		cnt += 1
		# Break out of the loop if too many iterations.
		if(cnt > 500)
			println("\nExceeded max iterations in run with seed ", seed, "\n")
			pass = false; break
		end
	end

	# Output to the data file as long as the simulation did not stall.
	if makeplots plotcircs(circvec, cnt, seed) end
	if pass
		@assert areafrac - get_afrac(circvec) < 100*eps(areafrac)
		plotcircs(circvec, -1, seed)
		datafile = string(geosfolder(), lpad(nbods,2,"0"), "-", seed, ".jld2")
		jldsave(datafile; circvec, areafrac, seed)
	end
end


# Call the main routine with given nbods and areafrac for 9 different seeds.
function make9geos(nbods::Int, areafrac::Float64)
	for ii=1:9
		make_geos(nbods, areafrac, ii, makeplots=false)
	end
end

end