# MAIN GOAL: Convert the lists of output text files to a Julia data file.
# Convention: nn indexes the timestep, el = 1:nbods indexes the bodies.

#--------------- BASIC STUFF ---------------#
using JLD2

# Later, once Jake makes the main Erosion module, change this to use that module.
include("erosion.jl")
# using

data_set() = "../output_data/erode-a06/"
savefile(label::AbstractString) = string(data_set(),"data-",label,".jld2")
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

# getparams: Get the parameters from a params file.
function getparams(paramvec::Vector, infovec::Vector)
	# Read from paramvec.
	circfile = String(paramvec[1])
	npts = paramvec[2]
	ibary, ifmm, ibc = Int(paramvec[3]), Int(paramvec[4]), Int(paramvec[5])
	epsfac, sigfac, dt, dtout = paramvec[6:9]
	fixpdrop, fixarea = Bool(paramvec[10]), Bool(paramvec[11])
	tfin = paramvec[12]
	maxl, nouter = Int(paramvec[13]), Int(paramvec[14])
	# Read from infovec.
	cntout = infovec[1]
	lastfile = Int(infovec[2])

	# Save the parameters in the updated object ParamSet.
	#= Note: due to the change in the file/folder labeling, 
	the infolder and label are not quite right, but these can be set manually. =#
	params = ParamSet(infolder=circfile, label="NA",
				npts=npts, ibary=ibary, ifmm=ifmm, ibc=ibc, 
				epsfac=epsfac, sigfac=sigfac, dt=dt, outstride=cntout,
				fixpdrop=fixpdrop, fixarea=fixarea, tfin=tfin,
				maxl=maxl, nouter=nouter)
	return params, lastfile
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
	# Get the basic meta-data.	
	paramvec = readvec( string(datafolder, "aparams.txt"))
	infovec = readvec( string(datafolder, "apinfo.txt"))
	params, lastfile = getparams(paramvec, infovec)

	println(params)

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