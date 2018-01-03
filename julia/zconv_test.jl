# Test the 2nd order convergence in dt.
# In typically use 01circ256, with tfin = 1e-2 and nits = 4 or higher.
# Vanishing time for 01circ1024aa is 0.017938. So 1e-2 is 56% of vanishing time.
# The runtime for nits = 4 is about 40 seconds.

include("main.jl")

#= Calculate the order of convergence given a vector of errors. =#
function order(Nfac::Integer, errv::Vector{Float64})
	return log(errv[1:end-1]./errv[2:end])./log(Nfac)
end
# Test the convergence wrt dt.
function convtest(tfin::Float64, nits::Int)
	paramsfile = "params.in"
	# Initialize arrays.
	dragv = zeros(Float64,nits)
	dtv = zeros(Float64,nits)
	ttv = zeros(Float64,nits)
	# Loop over different dt values.
	dt = 2e-3
	for nn=1:nits
		println("\n\nCONVERGENCE TEST: running number ", nn, " out of ", nits)
		thld_old = thlenden ???
		
		thlenden,params,tt = erosion(paramsfile,dt,tfin)
		dragv[nn] = drag(thlenden,params)[1]

		#
		dtv[nn] = dt
		ttv[nn] = tt
		dt = 0.5*dt
	end
	# Compute the error and order of convergence.
	errv = abs(dragv[1:end-1] - dragv[2:end])
	orderv = order(2,errv)
	# Print output to the screen.
	println("\n\nCONVERGENCE TEST RESULTS")
	println("\nThe final times: ", round(ttv,6))
	println("The dt values: ", round(dtv,6))
	println("\nThe drag values: ", round(dragv,3))
	println("The error: ", round(errv,4))
	println("The order: ", round(orderv,2))
end
