# If the folder exists, delete it. Then create a new folder.
function newfolder(foldername::AbstractString)
	if isdir(foldername)
		rm(foldername; recursive=true)
	end
	mkdir(foldername)
	return
end
# Plot the circles.
function plotcircs(circvec::Vector{CircType}, 
		figfolder::AbstractString, nfile::Int)
	# Parameters.
	npts = 128
	nfilestr = lpad("nfile",4,"0")
	figname = string(figfolder,"/circ", nfilestr,".pdf")
	width = 400; height = 400
	# Make the figure.
	pp = plot(xlim=(-1,1), ylim=(-1,1), size=(width,height), leg=false)
	# Plot each body.
	nbods = length(circvec)
	alpha = getalpha(npts)
	for nn = 1:nbods
		circ = circvec[nn]
		xx = circ.xc .+ circ.rad*cos.(alpha)
		yy = circ.yc .+ circ.rad*sin.(alpha)
		plot!(pp,xx,yy,color="black")
	end
	# Save the figure in a file.
	savefig(pp, figname)
	return
end

# Save the circle data to a file: radii and centers.
function save_circ_data(circvec::Vector{CircType})
	nbods = length(circvec)
	circfile = string(geosfolder(),nbods,"circ.circ")
	circdata = zeros(Float64, 0)
	radvec = zeros(Float64, nbods)
	for nn=1:nbods
		circ = circvec[nn]
		append!(circdata, [circ.rad, circ.xc, circ.yc])
		radvec[nn] = circ.rad
	end
	areafrac = 0.25*pi*sum(radvec.^2)
	data = [string("# nbods = ", nbods, ", areafrac = ", areafrac); nbods;
		"# data below: radius, xc, yc for each body."; circdata]
	writedata(data, circfile)
end
# Return thlen data for a circle of given radius and center.
function circthlen(npts::Int, rad::Float64, xc::Float64, yc::Float64)
	thlen = new_thlen()
	alpha = getalpha(npts)
	# Get the tangent angle, theta, and the total arclength, len.
	thlen.theta = 0.5*pi .+ alpha
	thlen.len = 2*pi*rad
	# Save the surface mean values, xsm and ysm.
	thlen.xsm = xc; thlen.ysm = yc
	return thlen
end
# Convert the circle data to theta-L data.
function save_thlen(circfile::AbstractString, thlenfile::AbstractString, npts::Int)
	# Read the circle data.
	circdata = readvec(circfile)
	nbods = round(Int, circdata[1])
	deleteat!(circdata,1)
	# Create the data vector to save the thlen values.
	thlendata = zeros(Float64, 2)
	thlendata[1] = npts
	thlendata[2] = nbods
	for nn = 1:nbods
		rad, xc, yc = circdata[1], circdata[2], circdata[3]
		thlen = circthlen(npts, rad, xc, yc)
		append!(thlendata, [thlen.theta; thlen.len; thlen.xsm; thlen.ysm])
		deleteat!(circdata,1:3)
	end
	writedata(thlendata, thlenfile)
end
# Convert the circle data to theta-L data.
function save_thlen(name::AbstractString, npts::Int)
	circfile = string(geosfolder(),name,".circ")
	thlenfile = string(geosfolder(),name,npts,".thlen")
	save_thlen(circfile, thlenfile, npts)
end

function geosfolder()
	return "../input_geos/"
end