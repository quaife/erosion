# Prep the data in Julia to view in Veusz.
# THIS CODE PROBABLY NEEDS AN OVERHAUL
include("main.jl")

function prep_data(datafolder::AbstractString)
	# Open the basic files.
	datafolder *= "/"
	paramsfile = string(datafolder, "aparams")
	infofile = string(datafolder,"apinfo.txt")
	# Get the meta-data.
	infovec = readvec(infofile)
	lastfile = Int(infovec[2])
	params = getparams(paramsfile)
	# Loop through the files inside folder.
	tm, por, perm, anis = [[] for ii=1:4]
	for nn=0:lastfile
		# Get the number of points and bodies.
		thlenden,nstr = get_thlenden(datafolder,nn)
		npts,nbods = getnvals(thlenden.thlenvec)
		# Read the geom data file to get the time.
		geomfile = string(datafolder,"geom",nstr,".dat")
		geomdata = readvec(geomfile)
		push!(tm, geomdata[1])
		# Read the areas data file to get the porosity.
		areafile = string(datafolder,"areas",nstr,".dat")
		nbods > 0 ? areavec = readvec(areafile) : areavec = [0]
		totarea = sum(areavec)
		push!(por, (4-totarea)/4)
		# Read the resistivity data file to get permeability and anisotropy.
		resfile = string(datafolder,"resistivity",nstr,".dat")
		resdata = readvec(resfile)
		nbods = resdata[1]
		res = resdata[2]
		resrot = resdata[3]
		push!(perm, 1/resdata[2])
		push!(anis, resdata[3]/resdata[2])
		# If nbods is zero, exit the loop.
		if nbods == 0 break end
	end
	# Clean up and write the output.
	tm = tm./tm[end]
	perm[end], anis[end] = Inf, NaN
	writedlm(string(datafolder,"aplots.txt"), [tm por perm anis])
end

prep_data("circ_shrink_20-2")
