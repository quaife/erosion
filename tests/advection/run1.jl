using Erosion

params = ParamSet(
		# The data file for the initial geometry
		infolder = "input_geos/",
		label = "100circ6",
		# Varied computational parameters
		npts = 128,                # The number of points per body, default 128.
		ibary = 1,                        # Use barycentric (1) or not (0).
		ifmm = 1,                        # Use the FFM (1) or not (0).
		ibc = 0,                        # Use slip BCs (1) or no-slip (0)
		# Varied physical parameters
		epsfac = 15,        # Smoothing parameter for curvature driven flow.
		sigfac = 10,        # Smoothing parameter for the stress.
		dt = 1e-2,                # The time step.
		outstride = 1,                # The number of steps per output, default 4.
		# Fixed physical parameters
		fixpdrop = 1,                # Fix the pressure drop (1) or not (0)
		fixarea = 0,                # Keep the area fixed (1) or not (0)
		tfin = 1e-2,                # The final time.
		# Fixed computational parameters
		maxl = 5000,                # The maximum number of GMRES iterations.
		nouter = 1024,                # The number of points on the outer boundary, default 1024.
)
main(params)

adv_params = AdvectionParamSet(
	jldfile = "data-02-1.jld",
	t = 0.0:0.01:3.0,
	npts = 121,
	thlenden_noutput = 4,
	# c_ana(x, t) = isBody(x[1], x[2], thlenden.thlenvec) ? 0.0 : tanh(1-x[1])^2
	ic = x -> x[1] == -1.0 ? 1.0 : 0.0 # fixed concentration left boundary
)
advect(adv_params)



# 100 bod streamlines
adv_params = AdvectionParamSet(
	jldfile = "../tests/advection/data-100circ6.jld",
	t = 0.0:0.01:3.0,
	npts = 121,
	thlenden_noutput = 42,
	ic = x -> x[1] == -1.0 ? 1.0 : 0.0 # fixed concentration left boundary
)

function porosity(thlenden::ThLenDenType)
	sum = 0
	for bod in thlenden.thlenvec
		xx, yy = Erosion.getxy(bod)
		sum += GeometryBasics.area(GeometryBasics.Point.(xx,yy))
	end
	return (4 - sum) / 4.0
end

# ANDRII
for i in 1:11:77

	thlenden = thlendens[i]
	
	plt = plot(aspect_ratio=1, size=(800,600))
	for bod in thlenden.thlenvec
		xx, yy = Erosion.getxy(bod)
		plot!(plt, Shape(xx, yy), label=:none, color=:black, xlims=(-1,1), ylims=(-1,1))
	end
	
	savefig(plt, "shape_$i.png")

end

# plot streamlines
plot_streamlines(thlenden.thlenvec, flow, 0.01, filename="run1_streamlines.png")

# solve and plot advection

    # characteristic solver
    function RK4(f, x, t, dt)
        k1 = dt * f(x, t)
        k2 = dt * f(x + k1/2.0, t + dt/2.0)
        k3 = dt * f(x + k2/2.0, t + dt/2.0)
        k4 = dt * f(x + k3, t + dt)
        x += 1 / 6 * (k1 + 2*k2 + 2*k3 + k4)
        return x
    end



bodyGrid = setup_isBody(X, thlenden.thlenvec) # initialize isBody matrix
state = advect(flow, concentration, bodyGrid)
plot_advection(state, X, "run1_advection.gif")


# check if point within body
function setup_isBody(X, thlenvec)
	
	function checkGridpoint(x,y)
		# check each body
		isBody = false
		for bod in thlenvec
			bx,by = Erosion.getxy(bod)
			j = length(bx) # start with line segment from bdy point N to bdy point 1
			for i in 1:length(bx)
				# check line segment from bdy point i to bdy point i+1
				if ( ((by[i] > y) != (by[j] > y)) && (x < ((bx[j] - bx[i]) * (y - by[i]) / (by[j] - by[i]) + bx[i])))
					isBody = !isBody
				end
				j = i
			end
			if isBody
				return isBody
			end
		end
		return isBody
	end

	return checkGridpoint.(getindex.(X,1), getindex.(X,2))

end

# streamline plot
function plot_streamlines(thlenvec, grid, dt; filename=nothing)

	# draw circles
	function plot_bodies(plt, thlenvec::Vector{ThetaLenType})	
		for ii = 1:lastindex(thlenvec)
			thlen = thlenvec[ii]
			xx,yy = Erosion.getxy(thlen)
			Plots.plot!(plt, xx, yy, color="black")
		end
	end

	# get grid and spacing
	x = unique(grid.xtar); y = unique(grid.ytar)
	nx = length(x); ny = length(y)
	dx = 2/(nx-1); dy = 2/(ny-1)
	
	# reshape list of points as 2d matrix
	U = reshape(grid.utar, ny, nx)
	V = reshape(grid.vtar, ny, nx)

	# build bicubic interpolation function on velocity
	uitp = interpolate(U', BSpline(Cubic(Line(OnGrid()))))
	u_interp = scale(uitp, minimum(x):dx:maximum(x), minimum(y):dy:maximum(y))
	vitp = interpolate(V', BSpline(Cubic(Line(OnGrid()))))
	v_interp = scale(vitp, minimum(x):dx:maximum(x), minimum(y):dy:maximum(y))

	# compute streamlines
	streams = []
	for (i, yi) in enumerate(y)

		@info "Body $i"
		
		# starting points on left boundary
		xi = x[1]
		stream_x = zeros(0); stream_y = zeros(0)

		# solve w/ forward euler until tracer hits a boundary
		its = 0
		while(xi <= maximum(x) && minimum(y) <= yi <= maximum(y) && its < 1000)
			its += 1
			# could also hit a body if the dt is too big 
			# if isBody(xi, yi, thlenvec)
			# 	break
			# end

			# store previous position
			push!(stream_x, xi); append!(stream_y, yi)
			
			# forward euler step
			try
				xi += u_interp(xi, yi) * dt
				yi += v_interp(xi, yi) * dt
			catch
				break # hit a boundary
			end
		
		end

		push!(streams, (stream_x, stream_y))

	end

	# plot
	plt = Plots.plot(size=(1800,1800), xlim=(minimum(x),maximum(x)), ylim=(minimum(y),maximum(y)), leg=false, aspect_ratio=:equal)
	for stream in streams
		Plots.plot!(plt, stream[1], stream[2], color="blue")
	end
	plot_bodies(plt, thlenvec)
	
	# show or save
	if filename === nothing
		display(plt)
	else
		Plots.savefig(plt, filename)
	end

end

# advection solver
function advect(flow::TargetsType, concentration, thlenvec::Vector{ThetaLenType})

	anim = @animate for (i, ti) in enumerate(t[1:end-1])

        # build bicubic interpolation function to find concentration off-grid
		# citp = interpolate(concentration, BSpline(Linear()))
		citp = interpolate(concentration, BSpline(Cubic(Reflect(OnGrid()))))
		c_interp = scale(citp, minimum(x):dx:maximum(x), minimum(y):dy:maximum(y))

		# characteristic equation dx/dt = -u
		dxdt(x, t) = [-u(x[1], x[2]), -v(x[1], x[2])]

		# solve characteristic equation
		sol = map(e -> solver(dxdt, e, ti, dt), X)

		# find characteristic tails outside domain or in solid bodies
		in_bounds = findall(e -> e[1] >= -1 && e[1] <= 1 && e[2] >= -1 && e[2] <= 1 && !isBody(e[1], e[2], thlenvec), sol)
		outside_domain = findall(e -> e[1] < -1 || e[1] > 1 || e[2] < -1 || e[2] > 1, sol)
		head_inside_bodies = findall(e -> isBody(e[1], e[2], thlenvec), X)
		tail_inside_bodies = findall(e -> isBody(e[1], e[2], thlenvec), sol)

		# upate concentration at gridpoints - set equal to conc at tail of characteristic
		concentration[in_bounds] = map(e -> c_interp(e[1], e[2]), sol[in_bounds])

		# don't advect if tail ends in body
		# concentration[tail_inside_bodies] .= itself, doy

		# don't advect at points inside bodies
		concentration[head_inside_bodies] .= 0

		# apply concentration boundary condition
		# fixed
		concentration[outside_domain] .= 1
		# or zero
		# concentration[outside_domain] .= 0
		# or decreasing over time
		# concentration[outside_domain] .= exp(-ti)
		# or periodic
		# concentration[outside_domain] .= abs(sin(pi/2*ti))

		# store state at this timestep
		push!(state, copy(concentration))

	end
	
	return state

end

# animated plot of advection, as heatmap
function plot_advection(state, X, plotfile)
	x = getindex.(X[:,1], 1)
	y = getindex.(X[1,:], 2)
	c_max = maximum(maximum.(state))
	c_min = minimum(minimum.(state))
	anim = @animate for concentration in state
		# heatmap(x, y, concentration, transpose=true, aspect_ratio=1, xlims=(-1,1), ylims=(-1,1), clims=(c_min, c_max))
		heatmap(x, y, concentration, transpose=true, aspect_ratio=1, xlims=(-1,1), ylims=(-1,1))
	end
	gif(anim, plotfile, fps=15)
end