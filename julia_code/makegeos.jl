#= Codes to make the initial geometries. =#
# Example commands to run
# makegeos(10, 0.5, 4)
# To make multiple geometries run: make9geos(10, 0.5)


#--------------- INITIALIZATION ---------------#
using Distributions
using Random
using Plots
using LinearAlgebra
using DelimitedFiles
include("./basic.jl")
include("./thetalen.jl")
include("./ioroutines.jl")
mutable struct CircType
	rad::Float64; xc::Float64; yc::Float64
end
function geosfolder()
	return "../input_geos/"
end
function figfolder()
	return "../zFigsGeos/"
end
include("iogeos.jl")


#--------------- MAIN ROUTINE ---------------#
# Call the main routine with given nbods and areafrac for 9 different seeds.
function make9geos(nbods::Int,areafrac::Float64)
	for ii=1:9
		makegeos(nbods,areafrac,ii)
	end
end
# Main routine to make the geometry.
function makegeos(nbods::Int, areafrac::Float64, seed::Int=1)
	# Parameters.
	buff = 0.08
	bolap = 0.03
	pow = 0.5
	dt = 5e-2
	fthresh = 1e-8
	# Check that the desired area fraction is not too high.
	# 0.91 is the absolute upper bound. 
	@assert areafrac < 0.71
	# Seed the random number generator.
	Random.seed!(seed)
	# Create the list of random radii.
	dray = Rayleigh()
	dchi = Chi(4)
	drad = dchi
	radvec = rand(drad, nbods)
	# Rescale the radii to achieve desired area fraction.
	radvec *= sqrt( 4*areafrac/ (pi*sum(radvec.^2)) )
	# Chose the provisional centers.
	duni = Uniform(-1,1)
	xc = rand(duni, nbods)
	yc = rand(duni, nbods)
	# Create the list of circles
	circvec = [CircType(radvec[nn],xc[nn],yc[nn]) for nn=1:nbods]
	# Shift the centers until no overlap.
	newfolder(figfolder())
	cnt = 0
	fx,fy,foverlap = forcesum(circvec,pow,buff,bolap)
	while((cnt < 50 || foverlap > 1e-10) && cnt < 800)
		# Plot the circles.
		println("count = ", cnt)
		#plotcircs(circvec,cnt,seed)
		# Shift the circles.
		sigma = 0.15*exp(-0.5*cnt*dt)
		shiftcircs(circvec,fx,fy,dt,sigma)
		cnt += 1
		fx,fy,foverlap = forcesum(circvec,pow,buff,bolap)
		println("foverlap = ",foverlap)
	end
	#plotcircs(circvec,cnt,seed)
	plotcircs(circvec,-1,seed)
	# Output to data files.
	save_circ_data(circvec,seed)
	return
end

#--------------- SUPPORTING ROUTINES ---------------#
# The repulsive force on circ1 due to circ2.
function fcircs(circ1::CircType, circ2::CircType, pow::Float64, buff::Float64, )
	# The critical distance.
	dcrit = (1+buff)*(circ1.rad + circ2.rad)
	# The vector from circ2 center to circ1 center, and the distance.
	v21 = [circ1.xc - circ2.xc, circ1.yc - circ2.yc]
	dist = norm(v21)
	# Compute the magnitude of the repulsive force.
	fmag = 0.
	if dist < dcrit
		fmag = (1-dist/dcrit)^pow
	end
	# Return the repulsive force
	if dist < 100*eps(dcrit)
		return 0*v21
	else
		return fmag*v21/dist
	end
end
# The repulisve force due to the walls
function fwall(circ::CircType, pow::Float64, buff::Float64)
	rcrit = (1+4*buff)*circ.rad
	# Get the force component in the x or y direction.
	function fcomp(xx::Float64)
		ff = 0.
		if xx < -1+rcrit
			ff = 1 - (xx+1)/rcrit
		elseif xx > 1-rcrit
			ff = -1 + (1-xx)/rcrit
		end
		return abs(ff)^pow * sign(ff)
	end
	# Get the force component in each direction.
	fx = fcomp(circ.xc)
	fy = fcomp(circ.yc)
	return [fx, fy]
end
# The sum of all forces on each circle.
function forcesum(circvec::Vector{CircType}, pow::Float64, 
		buff::Float64, bolap::Float64)
	nbods = length(circvec)
	#ftot = [zeros(Float64,2) for nn=1:nbods]
	fx,fy,fxo,fyo = [zeros(Float64,nbods) for nn=1:4]
	# Compute the total forces.
	for nn = 1:nbods
		finc = [0.,0.]
		finco = [0.,0.]
		# Compute the forces due to the other circles.
		for mm = 1:nbods
			if mm != nn
				finc += fcircs(circvec[nn],circvec[mm],pow,buff)
				finco += fcircs(circvec[nn],circvec[mm],pow,bolap)
			end
		end
		# Compute the forces due to the walls.
		finc += fwall(circvec[nn],pow,buff)
		finco += fwall(circvec[nn],pow,bolap)
		fx[nn] += finc[1]
		fy[nn] += finc[2]
		fxo[nn] += finco[1]
		fyo[nn] += finco[2]
	end
	foverlap = max(norm(fxo,Inf),norm(fyo,Inf))
	return fx,fy,foverlap
end
# Shift the circles with given force.
function shiftcircs(circvec::Vector{CircType}, 
		fxv::Vector{Float64}, fyv::Vector{Float64}, 
		dt::Float64, sigma::Float64)
	nbods = length(circvec)
	rvec = randn(2*nbods)
	sdt12 = sigma*sqrt(dt)
	for nn = 1:nbods
		circvec[nn].xc += dt*fxv[nn] + sdt12*rvec[2*nn-1]
		circvec[nn].yc += dt*fyv[nn] + sdt12*rvec[2*nn]
	end
	return
end
