# Test the 2nd order convergence in dt.
# I typically use 01circ256, with nits = 3 or higher, dt0 = 1e-4, tfin = 1e-2.
# Vanishing time for 01circ1024aa is 0.017938. So 1e-2 is 56% of vanishing time.
# The runtime for nits = 3 is about 16 minutes.
# For the test in the paper, I used nits = 7.

include("main.jl")

#= Calculate the order of convergence given a vector of errors. =#
function order(Nfac::Real, errv::Vector{Float64})
	return log(errv[1:end-1]./errv[2:end])./log(Nfac)
end
#= Calculate the error in two shapes. =#
function shape_error(thld1::ThLenDenType, thld0::ThLenDenType)
	npts1,nbods1,xv1,yv1 = getnxy(thld1)
	npts0,nbods0,xv0,yv0 = getnxy(thld0)
	assert(npts1==npts0)
	err = 1./sqrt(npts1) * norm([xv1;yv1]-[xv0;yv0])
	return err
end

# Test the convergence wrt dt.
function convtest(nits::Int=4, dt0::Float64=2e-4, tfin::Float64=1e-2; 
		paramsfile::AbstractString="params_conv.in")
	# Initialize simulation variables.
	thld1,params = startup(paramsfile)
	dtrat = 2.
	# Initialize output arrays.
	dragv,dtv,tv,cputime = [zeros(Float64,nits) for ii=1:4]
	err_shape,err_drag = [zeros(Float64,nits-1) for ii=1:2]
	# Loop over different dt values.
	dt = dt0
	for nn=1:nits
		println("\n\nCONVERGENCE TEST: running number ", nn, " out of ", nits)
		thld0 = deepcopy(thld1)
		thld1,params,tv[nn],cputime[nn] = erosion(paramsfile,dt,tfin)
		dragv[nn] = drag(thld1,params)[1]
		# Compute various errors.
		if nn>1
			err_shape[nn-1] = shape_error(thld1,thld0)
			err_drag[nn-1] = abs(dragv[nn] - dragv[nn-1])
		end
		# Save the dt and time values, then increment dt.
		dtv[nn] = dt
		dt *= 1./dtrat
	end
	# Compute the drag error and order of convergence.
	order_drag = order(dtrat, err_drag)
	order_shape = order(dtrat, err_shape)
	# Print output to the screen.
	println("\n\nCONVERGENCE TEST RESULTS")
	println("\nThe final times: ", round(tv,6))
	println("The dt values: ", round(dtv,6))
	# Print drag error.
	println("\nThe drag values: ", round(dragv,3))
	println("The drag errors: ", round(err_drag,4))
	println("The drag order: ", round(order_drag,2))
	# Print shape error.
	println("\nThe shape errors: ", round(err_shape,4))
	println("The shape order: ", round(order_shape,2))
	println()
	return dtv,err_shape,order_shape,cputime
end