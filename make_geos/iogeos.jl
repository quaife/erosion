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
	nfilestr = lpad(nfile,4,0)
	figname = string(figfolder,"/circ", nfilestr,".pdf")
	# Make the figure.
	pp = plot()
	xlim(-1,1); ylim(-1,1)
	# Plot each body.
	nbods = endof(circvec)
	alpha = getalpha(npts)
	for nn = 1:nbods
		circ = circvec[nn]
		xx = circ.xc + circ.rad*cos(alpha)
		yy = circ.yc + circ.rad*sin(alpha)
		pp = oplot(xx,yy,"-")
	end
	# Save the figure in a file.
	savefig(pp, figname, width=400, height=400)
	return
end

# Save the circle data to a file: radii and centers.
function save_circ_data(circvec::Vector{CircType})
	nbods = endof(circvec)
	circfile = string("../geos2/",nbods,"circ.dat")
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

# Create a circle of given radius and center.
function circthlen(npts::Integer, rad::Float64, xc::Float64, yc::Float64)
	thlen = new_thlen()
	alpha = getalpha(npts)
	# Get the tangent angle, theta, and the total arclength, len.
	thlen.theta = 0.5*pi + 2*pi*alpha
	thlen.len = 2*pi*rad
	# Save the surface mean values, xsm and ysm.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end

# Convert the circle data to theta-L data.
function save_thlen(circfile::AbstractString, npts::Float64)
	circdata = readvec(circfile)
	nbods = round(Int, circdata[1])
	deleteat!(circdata,1)
	#for nn = 1:nbods
end
