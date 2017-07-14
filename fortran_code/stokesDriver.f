      program stokesDriver
      implicit real*8 (a-h,o-z)

c      parameter (ninner = 1024)
      parameter (ninner = 128)
      parameter (nbodies = 2)
      parameter (nouter = 2**12)
      parameter (ntargets = 40)
      parameter (maxbodies = 10)

      dimension centerx(maxbodies),centery(maxbodies)
      dimension radius(maxbodies),phi(maxbodies)
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension den(2*ninner*nbodies + 3*nbodies + 2*nouter)
      dimension shear_stress(ninner*nbodies)
      dimension pressure(ninner*nbodies)
      dimension drag(2*nbodies)

      dimension xtar(ntargets*ntargets)
      dimension ytar(ntargets*ntargets)
      dimension utar(ntargets*ntargets)
      dimension vtar(ntargets*ntargets)
      dimension press_tar(ntargets*ntargets)

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(ninner)

      centerx(1) = -0.0d0
      centery(1) = 0.0d0
      centerx(2) = 0.8d0
      centery(2) = -0.1d0
      radius(1) = 2.d-1
      radius(2) = 2.d-1
      phi(1) = 0.d0
      phi(2) = 0.d0


      do j = 1,nbodies
        do k = 1,ninner
          theta = dble(k-1)*dtheta
          var_rad = radius(j)
          x((j-1)*ninner+k) = centerx(j) + var_rad*
     $          (dcos(phi(j))*dcos(theta) + dsin(phi(j))*dsin(theta))
          y((j-1)*ninner+k) = centery(j) + 2.d0*var_rad*
     $          (-dsin(phi(j))*dcos(theta) + dcos(phi(j))*dsin(theta))
        enddo
      enddo


      xmin = -1.d0
      xmax = 2.d0
      ymin = -5.d-1
      ymax = 5.d-1
      nx = ntargets
      ny = ntargets
      dx = (xmax - xmin)/dble(nx-1)
      dy = (ymax - ymin)/dble(ny-1)

      icount = 0
      do j = 1,nx
        do k = 1,ny 
          icount = icount + 1 
          xtar(icount) = xmin + dble(j-1)*dx
          ytar(icount) = ymin + dble(k-1)*dy
        enddo
      enddo

      ifmm = 1
      call stokesSolver(ninner,nbodies,nouter,ifmm,x,y,den)
c     pass in number of points and x and y coordinates and return the
c     density function on the boundary

      call computeShearStress(ninner,nbodies,nouter,x,y,den,
     $    shear_stress)
c     pass in the density function and return the shear_stress

      call computePressure(ninner,nbodies,nouter,x,y,den,pressure)
c     pass in the density function and return the pressure

      call computeDrag(ninner,nbodies,x,y,
     $      shear_stress,pressure,drag)
c     pass in the shear_stress and pressure and return the drag

      call computeVelocityPressureTargets(ninner,nbodies,nouter,
     $    x,y,den,ntargets*ntargets,xtar,ytar,utar,vtar,press_tar)


      open(unit=1,file='output/den.dat')
      open(unit=2,file='output/shear_stress.dat')
      open(unit=3,file='output/pressure.dat')
      open(unit=4,file='output/drag.dat')
      open(unit=8,file='output/x.dat')
      open(unit=9,file='output/y.dat')
      do k = 1,2*ninner*nbodies + 2*nouter + 3*nbodies
        write(1,1000) den(k)
      enddo
      do k = 1,ninner*nbodies
        write(2,1000) shear_stress(k)
        write(3,1000) pressure(k)
        write(8,1000) x(k)
        write(9,1000) y(k)
      enddo
      do k = 1,2*nbodies
        write(4,1000) drag(k)
      enddo
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=8)
      close(unit=9)

      open(unit=1,file='output/xtar.dat')
      open(unit=2,file='output/ytar.dat')
      open(unit=3,file='output/utar.dat')
      open(unit=4,file='output/vtar.dat')
      open(unit=8,file='output/press_tar.dat')
      do k = 1,ntargets**2
        write(1,1000) xtar(k)
        write(2,1000) ytar(k)
        write(3,1000) utar(k)
        write(4,1000) vtar(k)
        write(8,1000) press_tar(k)
      enddo


 1000 format(E25.16)



      end

