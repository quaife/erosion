# OBJECTIVE: Convert the lists of output text files to a Julia data file.
# Convention: nn indexes the timestep; bod = 1:nbods indexes the bodies.

# TO CALL give instructions


#--------------- BASIC STUFF ---------------#
using Erosion
using Erosion.ThetaLen # Permits use of new_thlenden.
using Erosion.MakeGeos # Permits use of CircType and plotcircs.
using JLD2
using DelimitedFiles

data_set() = "output_data/afrac06-text-data/"
savefile(label::AbstractString) = string("output_data/raw_data-",label,".jld2")
#-------------------------------------------------#

#--------------- LITTLE IO ROUTINES ---------------#
#= Many of these routines were initially in other files, either basic.jl, main.jl, or ioroutines.jl.
They became obsolete in the main part of the code but are still needed here. =#

# Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end

# Create new instances of each type.
new_thlen() = ThetaLenType([],0.,0.,0.,0.)
new_thlenvec(nbods::Int) = [new_thlen() for bod=1:nbods]

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
	for bod = 1:nbods
		# Read theta, len, xsm, ysm.
		thlenvec[bod].theta = invec[1:npts]
		thlenvec[bod].len = invec[npts+1]
		thlenvec[bod].xsm = invec[npts+2]
		thlenvec[bod].ysm = invec[npts+3]
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

# getparams: Get the parameters from a params file.
function getparams(datafolder::AbstractString, datalabel::AbstractString)
	# Get the basic meta-data.	
	paramvec = readvec( string(datafolder, "aparams.txt"))
	infovec = readvec( string(datafolder, "apinfo.txt"))
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
	cpu_hours = Int(infovec[3])

	# Save the parameters in the updated object ParamSet.
	#= Note: due to the change in the file/folder labeling, 
	the infolder and label are not quite right, but these can be set manually. =#
	params = ParamSet(infolder=circfile, label=datalabel,
				npts=npts, ibary=ibary, ifmm=ifmm, ibc=ibc, 
				epsfac=epsfac, sigfac=sigfac, dt=dt, outstride=cntout,
				fixpdrop=fixpdrop, fixarea=fixarea, tfin=tfin,
				maxl=maxl, nouter=nouter)
	return params, lastfile, cpu_hours
end
#-------------------------------------------------#

#--------------- REMAKE THE OUTPUT DATA ---------------#
# Convert the lists of output text files to a Julia data file.
function remake_output_data(datafolder::AbstractString, datalabel::AbstractString)
	params, lastfile, cpu_hours = getparams(datafolder, datalabel)
	println(params)

	# Loop through the files inside folder.
	thldvec = Vector{ThLenDenType}(undef, lastfile+1)
	for nn = 0:lastfile
		println("nn = ", nn)
		thlenden = get_thlenden(datafolder, nn)
		thldvec[nn+1] = thlenden
	end
	# Save thldvec in a Julia data file.
	jldsave(savefile(datalabel); params, thldvec, cpu_hours)
end

# Dispatch to run on simply a label using the folder given by data_set.
function remake_output_data(label::AbstractString)
	datafolder = string(data_set(), label, "/")
	remake_output_data(datafolder, label)
end

#--------------- REMAKE THE INPUT GEOMETRY DATA ---------------#
# Convert the input geometry text files to Julia data file.
function remake_input_geos(afrac_folder::AbstractString, nbods::Integer, seed::Integer)
	infile = string("input_geos/", afrac_folder, "/", nbods, "-", seed, ".circ")
	invec = readvec(infile)
	nbods = popfirst!(invec)
	circvec = Vector{CircType}(undef, 0)
	for bod = 1:nbods
		rad = popfirst!(invec)
		xc = popfirst!(invec)
		yc = popfirst!(invec)
		push!(circvec, CircType(rad, xc, yc))
	end
	areafrac = get_afrac(circvec)
	plotcircs(circvec, -1, seed)
	datafile = string("input_geos/", lpad(nbods,2,"0"), "-", seed, ".jld2")
	jldsave(datafile; circvec, areafrac, seed)
end
