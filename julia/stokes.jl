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

# stokes_thl: Call the Stokes solver given the theta-len values.
function stokes_thl(npts::Integer, nbods::Integer, 
		thetas::Vector{Float64}, lens::Vector{Float64}, 
		xcs::Vector{Float64}, ycs::Vector{Float64})
	ntot = npts*nbods
	xv = zeros(Float64,ntot)
	yv = zeros(Float64,ntot)
	for nn=1:nbods
		n1 = npts*(nn-1)+1
		n2 = npts*nn
		xv[n1:n2],yv[n1:n2] = getxy(thetas[n1:n2],lens[nn],xcs[nn],ycs[nn])
	end
	tau = stokes(npts,nbods,xv,yv)
	return abs(tau)
end