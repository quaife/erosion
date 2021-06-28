#= OBJECTIVE: Compute the density function, the shear stress, 
and other quantities by calling the Fortran code. =#

# Convention: bod = 1:nbods indexes the bodies; mm indexes the target points.

# The time-stepping routine only calls getstress!
# The erosion simulation only calls compute_density!
# The post-processing routine calls computestress, compute_qoi_targets!, regulargrid.

module DensityStress
export compute_density!, getstress!, computestress, compute_qoi_targets!, regulargrid

using ..SpectralMethods: gaussfilter
using ..ThetaLen: ParamSet, ThLenDenType, getnxy



#-------------------------------------------------#

using Statistics: mean

# QUESTION: Should I 'expose' this from the module ThetaLen???
using Parameters: @unpack	

#------------ Routines to compute the density function ------------#
# Little routine to rotate the x and y coordinates by 90 degrees CCW.
xyrot(xv::Vector{Float64}, yv::Vector{Float64}) = -yv, xv

# Wrapper to call Fortran routine stokessolver that computes the density function.
function compute_density(xx::Vector{Float64}, yy::Vector{Float64}, nbods::Int, params::ParamSet)
	@unpack npts, nouter, ifmm, ibary, ibc, maxl = params
	density = zeros(Float64, 2*npts*nbods + 3*nbods + 2*nouter)
	nits = zeros(Int,1)
	# Call the Fortran routine StokesSolver.
	ccall((:stokessolver_, "libstokes.so"), Nothing, 
		(Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},
		Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int}), 
		npts, nbods, nouter, ifmm, ibary, ibc, maxl, xx, yy, density, nits)
	println("The total number of GMRES iterations is ", nits[1],"\n\n")
	return density
end

# Dispatch for thlenden that computes the density function and saves in thlenden.
# This is function called by getstress!
function compute_density!(thlenden::ThLenDenType, params::ParamSet; rotation::Bool=false)
	density = rotation ? thlenden.denrot : thlenden.density
	if length(density) > 0; return; end
	println("Computing the density function with rotation = ", rotation)
	nbods,xv,yv = getnxy(thlenden)
	xv,yv = rotation ? xyrot(xv,yv) : (xv,yv)
	density = compute_density(xv,yv,nbods,params)
	rotation ? thlenden.denrot = density : thlenden.density = density 
end
#-----------------------------------------------------------------------#


#--------------- Routines to keep the pressure-drop fixed ---------------#
#= Note: In order to eventually calculate the stress, we must evaluate pressure on
a grid in order to enforce the condition of fixed pressure drop. =#

#--------------- Targets data type ---------------#
# Data type taht includes x-y coordinates of target points and u,v,pressure.
# Used in wrappers and in postprocess.
mutable struct TargetsType
	xtar::Vector{Float64}; ytar::Vector{Float64};
	utar::Vector{Float64}; vtar::Vector{Float64};
	ptar::Vector{Float64}; vortar::Vector{Float64};
end

#----- The quantities of interest: velocity, pressure, vorticity -----#
# Fortran wrapper to compute the QOI.
function compute_qoi_targets(xx::Vector{Float64}, yy::Vector{Float64},
		density::Vector{Float64}, xtar::Vector{Float64}, ytar::Vector{Float64},
		npts::Int, nbods::Int, nouter::Int, ibary::Int)
	ntargets = length(xtar)
	utar,vtar,ptar,vortar = [zeros(Float64,ntargets) for ii=1:4]
	ccall((:computeqoitargets_, "libstokes.so"), Nothing,
		(Ref{Int}, Ref{Int}, Ref{Int}, Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64},
		Ref{Int}, Ref{Float64}, Ref{Float64}, 
		Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
		npts, nbods, nouter, ibary, xx, yy, density, 
		ntargets, xtar, ytar, utar, vtar, ptar, vortar)
	return utar,vtar,ptar,vortar
end

# Dispatch for ThLenDenType and TargetsType to compute the QOI
function compute_qoi_targets!(thlenden::ThLenDenType, targets::TargetsType, params::ParamSet;
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods, xv, yv, density = getnxyden(thlenden, params, fixpdrop, rotation)
	targets.utar, targets.vtar, targets.ptar, targets.vortar = 
			compute_qoi_targets(xv,yv,density,targets.xtar,targets.ytar,npts,nbods,nouter,ibary)
	return
end

#--------------- The regular grid ---------------#
# Set up target points on a regular grid; return x and y.
function regulargrid(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	mx = length(xlocs)
	my = length(ylocs)
	ntargs = mx*my
	xtar = zeros(Float64, ntargs)
	ytar = zeros(Float64, ntargs)
	for mm = 1:mx
		m1 = my*(mm-1)+1
		m2 = my*mm
		xtar[m1:m2] .= xlocs[mm]
		ytar[m1:m2] .= ylocs
	end
	return xtar, ytar
end

# Set up target points on a regular grid; return targets.
function regulargridtargs(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	xtar,ytar = regulargrid(xlocs,ylocs)
	targets = TargetsType([], [], [], [], [], [])
	targets.xtar = xtar
	targets.ytar = ytar
	return targets
end

#--------------- Enforce the pressure drop condition ---------------#
#= Calculate the pressure drop from -x0 to x0. Also get the average flux while at it. 
Note: this routine assumes that umax = 1; If different, need to apply rescaling. =#
function getpdrop(thlenden::ThLenDenType, params::ParamSet, x0::Float64 = 2.0; rotation::Bool=false)
	# Set up targets points on two vertical slices.
	nypts = 13
	dy = 2/nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = regulargridtargs([x0],ylocs)
	tarm = regulargridtargs([-x0],ylocs)
	compute_qoi_targets!(thlenden, tarp, params, fixpdrop=false, rotation=rotation)
	compute_qoi_targets!(thlenden, tarm, params, fixpdrop=false, rotation=rotation)
	# Compute the pressure drop.
	pplus = mean(tarp.ptar)
	pminus = mean(tarm.ptar)
	pdrop = pminus-pplus
	# Compute the flux.
	qplus = mean(tarp.utar)
	qminus = mean(tarm.utar)
	qavg = 0.5*(qplus+qminus)
	#= The discharge should be exactly the same at any location x.
	So check that it is the same at x0 and -x0. =#
	qreldiff = (qplus-qminus)/qavg
	if qreldiff > 1e-3
		@warn string("The flux does not match at x0 and -x0: qreldiff = ", round(qreldiff,sigdigits=3))
	end
	return pdrop, qavg
end

# Get umax to rescale the density function.
function getumax(thlenden::ThLenDenType, params::ParamSet, fixpdrop::Bool, x0::Float64 = 2.0)
	# NOTE: With u = 1-y^2 and x0 = 2, the pressure drop is pdrop = 8.
	# This code should work for pipe flow or uniform flow: ibc = 0 or 1.
	umax = 1.
	if fixpdrop
		pdrop = getpdrop(thlenden, params)[1]
		umax = 4*x0/(pdrop + 1e-4)
		println("Fixing pdrop, umax = ", round(umax,sigdigits=3))
	end
	return umax
end

# Get the density function, and a few other things, depending on fixpdrop and rotation.
# This function is called by compute_stress. 
#= Note: The density should have already been computed. 
This routine simply accesses it and rescales if appropriate. =#
function getnxyden(thlenden::ThLenDenType, params::ParamSet, fixpdrop::Bool, rotation::Bool)
	nbods, xv, yv = getnxy(thlenden)
	# Consider fixpdrop.
	rescale = getumax(thlenden, params, fixpdrop)
	# Consider rotation.
	if rotation
		xv,yv = xyrot(xv,yv)
		density = rescale * thlenden.denrot
	else
		density = rescale * thlenden.density
	end
	return nbods, xv, yv, density
end
#-----------------------------------------------------------------------#


#------------ Routines to compute the shear stress ------------#
# Fortran wrapper that computes the shear stress.
function compute_stress(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int, ibary::Int)
	tau = zeros(Float64, npts*nbods)
	if nbods > 0
		ccall((:computeshearstress_, "libstokes.so"), Nothing,
			(Ref{Int},Ref{Int},Ref{Int},
			Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int},Ref{Float64}),
			npts, nbods, nouter, xx, yy, density, ibary, tau)
	end
	return tau
end

# Dispatch for ThLenDenType to compute the shear stress.
function compute_stress(thlenden::ThLenDenType, params::ParamSet; 
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods, xv, yv, density = getnxyden(thlenden,params,fixpdrop,rotation)
	tau = compute_stress(xv,yv,density,npts,nbods,nouter,ibary)
	return reshape(tau, npts, nbods)
end

#Compute the smoothed absolute stress atau and saves it in thlenden.thlenvec.atau.
# This is the main routine called by the erosion simulation in timestepping.
function getstress!(thlenden::ThLenDenType, params::ParamSet)
	den_time = @elapsed		compute_density!(thlenden, params)
	str_time = @elapsed		stress = compute_stress(thlenden, params, fixpdrop=params.fixpdrop)
	
	println("\nTime taken to compute density = ", round(den_time, sigdigits=3), "sec.")
	println("Time taken to compute stress = ", round(str_time, sigdigits=3), "sec.")
	
	# Smooth the stress and also save the mean of smoothed stress atau.
	for bod = 1:length(thlenden.thlenvec)
		stress[:,bod] = gaussfilter( abs.(stress[:,bod]), params.sigma)
		thlenden.thlenvec[bod].matau = mean(stress[:,bod])
	end
	return stress
end

end