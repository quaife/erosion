# postporcess.jl
# Post-processing routines.

# postprocess: Use the saved data to compute stuff.
function postprocess(foldername::AbstractString)
    # Define the data folder.
    datafolder = string("../datafiles/",foldername,"/")
    # Read the params data file.
    paramsfile = string(datafolder,"params.dat")
    paramvec = readvec(paramsfile)
    nouter = paramvec[2]
    ntimes = paramvec[end]
    # Read the data at each time step.
    for cnt=0:ntimes
        # Get the file name at each time.
        cntstr = lpad(cnt,4,0)
        geomfile = string(datafolder,"geom",cntstr,".dat")
        densityfile = string(datafolder,"density",cntstr,".dat")
        # Extract thlenvec and density.
        tt,thlenvec = read_geom_file(geomfile)
        density = readvec(densityfile)
        thlenden = new_thlenden(thlenvec,density)

        #--------------------------------------#
        # Compute velocity and pressure at a set of target points.
        xlocs = [-2.5, -2.0, -1.5, 1.5, 2.0, 2.5]
        ylocs = collect(-0.8: 0.2: 0.8)
        targets = setuptargets(xlocs,ylocs)
        compute_velpress_targets!(thlenden,targets,nouter)
        # Save the output to a data file.
        targfile = string(datafolder,"targs",cntstr,".dat")
        iostream = open(targfile, "w")
        label = string("# Stuff")
        writedlm(iostream, [label; targets.xtar; targets.ytar; 
            targets.utar; targets.vtar; targets.ptar])
        close(iostream)

        #--------------------------------------#
        # Compute the permeability
        rt1,rb1 = resistivity(thlenden, nouter, 1.5)
        rt2,rb2 = resistivity(thlenden, nouter, 2.0)
        rt3,rb3 = resistivity(thlenden, nouter, 2.5)

        #--------------------------------------#
        # Compute the drag on each body.
        # Note: the stress is not smoothed and absolute value is not taken.
        pressvec = compute_pressure(thlenden,nouter)
        tauvec = compute_stress(thlenden,nouter)
        npts,nbods = getnvals(thlenvec)
        dragxvec,dragyvec = [zeros(Float64,nbods) for ii=1:2]
        for nn=1:nbods
            # Get the pressure and stress and boddy nn.
            n1,n2 = n1n2(npts,nn)
            press = pressvec[n1:n2]
            tau = tauvec[n1:n2]
            # Get the tangent/normal vectors and arc length increment.
            sx,sy,nx,ny = getns(thlenvec[nn].theta)
            ds = thlenvec[nn].len / npts
            # Compute the drag force.
            # Note: I believe both should be plus signs due to the conventions.
            dragxvec[nn] = sum(press.*nx + tau.*sx)*ds
            dragyvec[nn] = sum(press.*ny + tau.*sy)*ds
        end
        # Save the output to a data file.
        dragfile = string(datafolder,"drag",cntstr,".dat")
        iostream = open(dragfile, "w")
        label = string("# Drag data for ",nbods," bodies. All xdrags then all ydrags.")
        writedlm(iostream, [label; dragxvec; dragyvec])
        close(iostream)
    end
    return
end

# setuptargets: Set up the target points.
function setuptargets(xlocs::Vector{Float64}, ylocs::Vector{Float64})
    nx = endof(xlocs)
    ny = endof(ylocs)
    ntargs = nx*ny
    xtar = ones(Float64,ntargs)
    ytar = ones(Float64,ntargs)
    for nn=1:nx
        n1,n2 = n1n2(ny,nn)
        xtar[n1:n2] = xlocs[nn]
        ytar[n1:n2] = ylocs
    end
    targets = TargetsType(evec(), evec(), evec(), evec(), evec())
    targets.xtar = xtar
    targets.ytar = ytar
    return targets
end

# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64)
    # Set up targets points on a y-grid for midpoint rule.
    nypts = 13
    dy = 2./nypts
    ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
    # Target points for plus/minus x0.
    tarp = setuptargets([x0],ylocs)
    tarm = setuptargets([-x0],ylocs)
    # Comopute the velocities and pressures on each set of target points.
    compute_velpress_targets!(thlenden,tarp,nouter)
    compute_velpress_targets!(thlenden,tarm,nouter)
    # Compute the cross-sectional average pressure and discharge
    pplus = mean(tarp.ptar)
    pminus = mean(tarm.ptar)
    qplus = mean(tarp.utar)
    qminus = mean(tarm.utar)
    qavg = 0.5*(qplus+qminus)
    #= The discharge should be exactly the same at any location x.
    So check that it is the same at x0 and -x0. =#
    qreldiff = (qplus-qminus)/qavg
    assert(qreldiff < 1e-6)
    # Calculate the total resistivity
    rtot = (pminus - pplus)/(2*x0*qavg)
    # Calculate the resisitvity due only to the bodies.
    rbods = x0*(rtot - 3)
    # For testing.
    #println("At x0 = ", x0, " the total resistivity is: ", signif(rtot,3))
    #println("At x0 = ", x0, " the matrix resistivity is: ", signif(rbods,3))
    return rtot,rbods
end
