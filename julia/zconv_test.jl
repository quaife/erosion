2# Test the 2nd order convergence in dt.
# In params.dat, I use 01circ256.in; tfin = 0.01.
include("main.jl")
include("postprocess.jl")

#= Calculate the order of convergence given a vector of errors. =#
function order(Nfac::Integer, errv::Vector{Float64})
	return log(errv[1:end-1]./errv[2:end])./log(Nfac)
end

# Test the convergence wrt dt.
function convtest()
	nits = 6
	dragv = zeros(Float64,nits)
	dtv = zeros(Float64,nits)
	ttv = zeros(Float64,nits)
	# Loop over different dt values.
	dt = 2e-3
	for nn=1:nits
		thlenden,params,tt = erosion(dt)
		dragv[nn] = drag(thlenden,params.nouter)[1]
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

convtest()
