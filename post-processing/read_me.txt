
(1) stokesDriver_geo_den_target.f:
    Use geom files and density files to find the interest parameters (velocity, vorticity, pressure etc.) at target point. 

(2) stokesDriver_geo_den_pressure_drop.f: 
    Use geom files and density files to find the pressures at inlet and outlet.

(3) stokesTracer.f: 
    Use geom files and density files to find the streamlines by RK4.
    
(4) stokesTracer_long.f: 
    Make the long streamlines for calculting anomalous dispersion rate. 
    
-------------------------------------------------------------------------------------------    

(5) vorticity_target.m: [do it after (1) & (2)]
    Do the pressure drop on velocity field and vorticity field, then make a vedio. 

(6) porosity.m: [do it after (1)]
    find the porosity of the geometric structure.

(7) tortuosity_Eulerian.m: [do it after (1)]

(8) tort_porosity_fitting.m: [do it after (1), (6), & (7)]
    find the fitting curves for tortuosity.

(9) tortuosity_Lagragian.m: [ do it after (3) ]

(10) tracerplot.m: [do it atfer (3)]
    tracer vedio.
    
(11) anomalous_diffusion_long.m: [do it after (4)]
     get the second moment (anomalous dispersion rate) of the streamlines.
     
(12) gap_triangulation.m: 
     Find the Delaunay Triangulation of our grain distribution.

(13) gap_draw.m: [for figure 3]
     draw a simple of Delaunay Triangulation and a corresponding distribution of pore size.

(14) hist_fitting.m: [do it after (12)]
     draw the distribution of pore size with a fitting Weibull distribution 
     
(15) shape_eroding.m: [using function color_line3.m]
     draw the shear stress on the sigle body. 
     
    
