#= stokes.jl: Functions involving the Stokes solver.
All routines work for multiple bodies. =#

# stokes: Julia wrapper to call the Fortran stokessolver
function stokes(npts::Integer, nbods::Integer, xx::Vector{Float64}, yy::Vector{Float64})
	ntot = npts*nbods
	tau = zeros(Float64, ntot)
	ccall((:stokessolver_, "libstokes.so"),
		Void, (Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
		&npts, &nbods, xx, yy, tau)
	return tau
end
#= stokes!: Dispatch for vector of ThetaLenType. Updates atau in each thlen.
Note: All of the entries in thlenv should already be loaded with xx and yy values. =#
function stokes!(thlenv::Vector{ThetaLenType})
	nbods = endof(thlenv)
	npts = endof(thlenv[1].theta)
	ntot = nbods*npts
	xv = zeros(Float64,ntot); yv = zeros(Float64,ntot)
	# Put all of the xy values in a single vector.
	for nn = 1:nbods
		# Compute the xy coordinates if they are not already loaded in thlen.
		getxy!(thlenv[nn])
		# Put the values into a single vector.
		n1 = npts*(nn-1)+1
		n2 = npts*nn
		xv[n1:n2], yv[n1:n2] = thlenv[nn].xx, thlenv[nn].yy
	end
	# Call the stokessolver.
	tau = stokes(npts,nbods,xv,yv)
	# Update the atau value in each of the thlen variables.
	for nn = 1:nbods
		n1 = npts*(nn-1)+1
		n2 = npts*nn
		thlenv[nn].atau = abs(tau[n1:n2])
	end
	return
end
