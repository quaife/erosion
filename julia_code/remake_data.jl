# MAIN GOAL: Convert the lists of output text files to a Julia data file.
# Convention: nn indexes the timestep, el = 1:nbods indexes the bodies.

#--------------- BASIC STUFF ---------------#
using JLD2

include("basic.jl")
data_set() = "../output_data/erode-a06/"
savefile(label::AbstractString) = string(data_set(),"data-",label,".jld2")

# Unused for now, but will be used if/when I delete ParamType
#include("run0.jl")
#-------------------------------------------------#





# basic.jl: Basic routines such as datatypes.
using DelimitedFiles
# readvec: Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end









#--------------- LITTLE IO ROUTINES ---------------#
#= These routines were initially in other files, either basic.jl, main.jl, or ioroutines.jl.
They became obselote in the main part of the code but are still needed 
to remake the data, so I moved them here.=#

# ParamType: The outdated datatype of parameters.
# SHOULD EVENTUALLY DELETE THIS DATATYPE AND USE THE NEW PARAMSET INSTEAD.
mutable struct ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; 
	nouter::Int; ifmm::Int; ibary::Int; ibc::Int;
	maxl::Int; fixarea::Bool; fixpdrop::Bool;
	npts::Int; tfin::Float64; cntout::Int; cput0::Float64;
	circfile::AbstractString; paramsfile::AbstractString
end

# getparams: Get the parameters from a params file.
function getparams(paramsfile::AbstractString)
	# Read the parameters.
	paramvec = readvec(string(paramsfile,".txt"))
	circfile = paramvec[1]
	npts = paramvec[2]
	ibary, ifmm, ibc = Int(paramvec[3]), Int(paramvec[4]), Int(paramvec[5])
	epsfac, sigfac, dt, dtout = paramvec[6:9]
	fixpdrop, fixarea = Bool(paramvec[10]), Bool(paramvec[11])
	tfin = paramvec[12]
	maxl, nouter = Int(paramvec[13]), Int(paramvec[14])
	# Calculate the needed quantities.
	epsilon = epsfac/npts
	sigma = sigfac/npts
	# Unused definitions.
	#cntout = max(round(Int,dtout/dt),1)
	#cput0 = time()
	cntout = 0; cput0 = 0.0

	# Save the parameters in an object.
	# SHOULD EVENTUALLY USE PARAMSET INSTEAD.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,ibary,ibc,maxl,
		fixarea,fixpdrop,npts,tfin,cntout,cput0,circfile,paramsfile)
	return params
end

# Create new instances of each type.
new_thlen() = ThetaLenType([],0.,0.,0.,0.)
new_thlenvec(nbods::Int) = [new_thlen() for el=1:nbods]

# read_geom_file: Read a geometry file.
#= The data in the file is t (physical time), npts, nbods, 
then theta, len, xsm, ysm, xx, yy for each body. =#
function read_geom_file(filename::AbstractString)
	invec = readvec(filename)
	tt = invec[1]
	npts = round(Int,invec[2])
	nbods = round(Int,invec[3])
	deleteat!(invec,1:3)
	@assert length(invec) == nbods*(3*npts+3)
	thlenvec = new_thlenvec(nbods)
	for el=1:nbods
		# Read theta, len, xsm, ysm.
		thlenvec[el].theta = invec[1:npts]
		thlenvec[el].len = invec[npts+1]
		thlenvec[el].xsm = invec[npts+2]
		thlenvec[el].ysm = invec[npts+3]
		deleteat!(invec,1:npts+3)
		# Reading xx and yy is now obselete, but I still need to delete entries in invec.
		deleteat!(invec,1:2*npts)
	end
	return tt, thlenvec
end

# read_density_file: Read the density file to get the density and rotated density.
function read_density_file(filename::AbstractString)
	dendata = readvec(filename)
	nn = length(dendata)
	@assert iseven(nn)
	nhalf = div(nn,2)
	density = dendata[1:nhalf]
	denrot = dendata[nhalf+1:nn]
	return density, denrot
end

# get_thlenden: Get thlenden from a datafolder.
function get_thlenden(datafolder::AbstractString, cnt::Int)
	# Get the file name at each time.
	cntstr = lpad(string(cnt),4,"0")
	geomfile = string(datafolder,"geom",cntstr,".dat")
	densityfile = string(datafolder,"density",cntstr,".dat")
	# Extract thlenvec, density, and denrot.
	tt, thlenvec = read_geom_file(geomfile)
	density, denrot = read_density_file(densityfile)
	# Create variable thlenden
	thlenden = new_thlenden(thlenvec)
	thlenden.tt = tt
	thlenden.density = density
	thlenden.denrot = denrot
	return thlenden
end
#-------------------------------------------------#


#--------------- REMAKE THE DATA ---------------#
# Convert the lists of output text files to a Julia data file.
function remake_data(datafolder::AbstractString, datalabel::AbstractString)
	# Open the basic files.
	paramsfile = string(datafolder, "aparams")
	infofile = string(datafolder,"apinfo.txt")
	# Get the meta-data.
	infovec = readvec(infofile)
	lastfile = Int(infovec[2])
	params = getparams(paramsfile)
	params.cntout = Int(infovec[1])
	params.cput0 = infovec[3]
	# Loop through the files inside folder.
	thldvec = Vector{ThLenDenType}(undef, lastfile+1)
	for nn = 0:lastfile
		println("nn = ", nn)
		thlenden = get_thlenden(datafolder, nn)
		thldvec[nn+1] = thlenden
	end
	# Save thldvec in a Julia data file.
	jldsave(savefile(datalabel); thldvec, params)
end

# Dispatch to run on simply a label using the folder given by data_set.
function remake_data(label::AbstractString)
	datafolder = string(data_set(), label, "/")
	remake_data(datafolder, label)
end
#-------------------------------------------------#

remake_data("20-2")