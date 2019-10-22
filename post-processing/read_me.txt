
(1) stokesDriver_geo_den_target.f:
    Use geom files and density files to find the interest parameters (velocity, vorticity, pressure etc.) at target point. 

(2) stokesDriver_geo_den_pressure_drop.f: 
    Use geom files and density files to find the pressures at inlet and outlet.

(3) stokesTracer.f: 
    Use geom files and density files to find the streamlines by RK4.

(4) vorticity_target.m: [do it after (1) & (2)]
    Do the pressure drop on velocity field and vorticity field, then make a vedio. 

(5) porosity.m: [do it after (1)]
    find the porosity of the geometric structure.

(6) tortuosity_Eulerian.m: [do it after (1)]

(7) tort_porosity_fitting.m: [do it after (1), (5), & (6)]
    find the fitting curves for tortuosity.

(8) tortuosity_Lagragian.m: [ do it after (3) ]

(9) tracerplot.m: [do it atfer (3)]
    tracer vedio.
