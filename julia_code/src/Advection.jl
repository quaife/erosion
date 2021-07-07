### NFG / OUT OF DATE ###
### THIS MODULE NEEDS TO BE UPDATED AFTER THE JULY 2021 MERGE ###

module Advection

    export advect, AdvectionParamSet

    # using JLD: load
    using Erosion.ThetaLen
    # using Erosion: main, compute_density, compute_density!, compute_qoi_targets!
    using Parameters: @with_kw, @unpack
    using Base.Iterators: product
    # using Plots
    using GeometricalPredicates: inpolygon, Polygon, Point
    using Interpolations
    import ScatteredInterpolation

    @with_kw struct AdvectionParamSet
        # required parameters
        jldfile::String
        t::Array{Float64}
        ic::Function
        npts::Int
        # defaulted parameters
        thlenden_noutput::Int=1
        # derived parameters
        dt::Float64 = t[2] - t[1]
        dx::Float64 = 2/(npts-1)
        dy::Float64 = 2/(npts-1)
    end

    function advect(adv_params::AdvectionParamSet)

        @unpack jldfile, t, ic, npts, thlenden_noutput, dt, dx, dy = adv_params

        # load erosion simulation
        @info "Loading erosion data..."
        jf = load(jldfile)
        params = jf["params"]
        thlenden = jf["thlendens"][thlenden_noutput]

        # set up regular grid
        xlocs = collect(-1:dx:1); ylocs = collect(-1:dy:1)
        xtar = [x for x in xlocs for y in ylocs]; ytar = [y for x in xlocs for y in ylocs]
        grid = TargetsType(xtar, ytar, [], [], [], [])
        
        # precompute density kernel
        @info "Finding density kernel..."
        @time compute_density!(thlenden, params)

        # precompute velocity on gridpoints
        @info "Precomputing on-grid velocity ..."
        @time compute_qoi_targets!(thlenden, grid, params)

        # precompute velocity at characteristic tails
        @info "Precomputing off-grid velocity..."
        @time tails = precompute_RK4(thlenden, params, grid, dt)

        # initial concentration
        c = reshape(mapslices(ic, hcat(grid.xtar, grid.ytar); dims=2), npts^2)


        ### TODO:
            # forget about the whole RBF interp on a whole body thing, we just have to fill out the few
            # gridpoints in a body's interior NEAR the very FEW tails that end so close to a body.
            # (1) above, get rid of excess tails that start and end inside bodies
            # (2) find tails within 2 gridpoints of a body boundary
            # (3) do a VERY local RBF extension with a few neighboring, exterior points
            # (4) be careful of case where two bodies very close, do not want to include these in interp
            # (5) track all of the local interpolations to do, so that re-interp and re-eval during time loop is v. fast
        
        # find tails that end near boundary
        @time begin
            tar = hcat(grid.xtar, grid.ytar)
            # get boundary
            bodies = [Polygon(Point.(xx,yy)...) for (xx,yy) in Erosion.getxy.(thlenden.thlenvec)]

        end

        # RBF extension into body






        ### CHECK IF TAILS NEED INTERP AT BOUNDARY

        ### IF THEY DO, THEN

        ### GET BODIES TO PROCESS
        
        bods = thlenden.thlenvec

        # local rbf interp around body
        for bod in bods
            @time begin
                tar = hcat(grid.xtar, grid.ytar)
                # get boundary
                xx, yy = Erosion.getxy(bod)
                poly = Polygon(Point.(xx,yy)...)
                # find reg gridpoints in box containing body
                box = findall(
                    (minimum(xx) - 2dx .<= tar[:,1] .<= maximum(xx) + 2dx) .&
                    (minimum(yy) - 2dy .<= tar[:,2] .<= maximum(yy) + 2dy)
                )
                # find reg gridpoints close to boundary (within 2 gridpoints of boundary)
                border = box[findall([minimum(hypot.(px .- xx, py .- yy)) <= 2dx for (px,py) in zip(tar[box,1], tar[box,2])])]
                # separate interior and exterior gridpoints
                ibd = findall(inpolygon.(Ref(poly), Point.(tar[border,1],  tar[border,2])))
                obd = setdiff(1:length(border[:,1]), ibd)
                ibd = border[ibd]; obd = border[obd]
                # interpolate on exterior points neighboring boundary
                @time itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), tar[obd,:]', c[obd])
                # fill in interior points neighboring boundary
                @time c[ibd] = ScatteredInterpolation.evaluate(itp, tar[ibd,:]')
            end

            # c = (x,y) -> inpolygon(poly, Point(x, y)) ? 0.0 : exp(-x^2) + exp(-y^2)
            # ic = [c(x,y) for (x,y) in zip(tar[box,1], tar[box,2])]
            # icobd = ic[map(x->x.start,searchsorted.(Ref(box), obd))]

            itp = ScatteredInterpolation.interpolate(ThinPlate(), tar[obd,:]', icobd)
            xtp = mapslices(p -> evaluate(itp, [p[1], p[2]]), tar[ibd,:], dims=2)
            ic[map(x->x.start,searchsorted.(Ref(box), ibd))] = xtp
            ### GET MATRIX SO WE CAN PRECOMPUTE THIS ONCE SINCE IT WILL BE TIME INDEPENDENT


            # heatmap(unique(grid.xtar[box]), unique(grid.ytar[box]), reshape(ic,(length(unique(grid.ytar[box])),length(unique(grid.xtar[box])))))

            # plot(Erosion.getxy(bod)...; label="", aspect_ratio=1)
            scatter!(box[:,1], box[:,2], label="", markersize=1)
            scatter!(border[:,1], border[:,2], label="")
            
        end

        # loop over time
        @info "Begin advection!"
        concentration = [c]
        @time for (i, ti) in enumerate(t[1:end-1])

            c_itp = interpolate((xlocs, ylocs), transpose(reshape(concentration[end], (npts,npts))), Gridded(Linear()))
            c_etp = extrapolate(c_itp, Flat())
            push!(concentration, c_etp.(tails[:,1], tails[:,2]))
            
        end
        
        # gif
        @info "Generating plots..."
        # @time anim = @animate for c in concentration
        #     heatmap(xlocs, ylocs, reshape(c, (npts,npts)); transpose=true, aspect_ratio=1, xlims=(-1,1), ylims=(-1,1))
        # end

        # display(gif(anim, fps=30))

        @info "Done!"

    end

    # time-independent RK4 solve
    function precompute_RK4(thlenden::ThLenDenType, params::ParamSet, k1_grid::TargetsType, dt::Float64)
        k2_grid = TargetsType(
            [x - u * dt / 2.0 for (x,u) in zip(k1_grid.xtar, k1_grid.utar)],
            [y - v * dt / 2.0 for (y,v) in zip(k1_grid.ytar, k1_grid.vtar)],
            [], [], [], []
        )
        compute_qoi_targets!(thlenden, k2_grid, params)
        k3_grid = TargetsType(
            [x - u * dt / 2.0 for (x,u) in zip(k1_grid.xtar, k2_grid.utar)],
            [y - v * dt / 2.0 for (y,v) in zip(k1_grid.ytar, k2_grid.vtar)],
            [], [], [], []
        )
        compute_qoi_targets!(thlenden, k3_grid, params)
        k4_grid = TargetsType(
            [x - u * dt for (x,u) in zip(k1_grid.xtar, k3_grid.utar)],
            [y - v * dt for (y,v) in zip(k1_grid.ytar, k3_grid.vtar)],
            [], [], [], []
        )
        compute_qoi_targets!(thlenden, k4_grid, params)

        utar = 1 / 6 * (k1_grid.utar .+ 2 * k2_grid.utar .+ 2 * k3_grid.utar .+ k4_grid.utar)
        vtar = 1 / 6 * (k1_grid.vtar .+ 2 * k2_grid.vtar .+ 2 * k3_grid.vtar .+ k4_grid.vtar)
        xtar = k1_grid.xtar .- dt * utar
        ytar = k1_grid.ytar .- dt * vtar

        return hcat(xtar, ytar)
    end

end