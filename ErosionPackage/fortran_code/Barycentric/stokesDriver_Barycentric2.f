      program stokesDriver
      implicit real*8 (a-h,o-z)

      parameter (ninner = 128)
      parameter (nbodies = 2)
      parameter (nouter = 2**8)
      parameter (ntargets = 120)
      parameter (maxbodies = 10)

c      parameter (ninnc = 2**8,nbeta = 1,ninner = nbeta*ninnc)      
c      parameter (noutc = 2**8, nouter = nbeta*noutc)      
      
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
      dimension vort_tar(ntargets*ntargets)
      
c     Testing for stress tensor       
      dimension xtar_test(ninner*nbodies)
      dimension ytar_test(ninner*nbodies)
      dimension utar_test(ninner*nbodies)
      dimension vtar_test(ninner*nbodies)
      dimension press_tar_test(ninner*nbodies)
      dimension vort_tar_test(ninner*nbodies)      
      
      
      dimension iside(ntargets*ntargets)
      dimension inear(ntargets*ntargets)

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      dtheta = twopi/dble(ninner)

!      centery(1) = 0.0d0
!      centerx(2) = 0.5d0
!      centerx(3) = -0.4d0
!      centery(2) = -0.1d0
!      centery(3) = 0.1d0
!      radius(1) = 2.d-1
!      radius(2) = 1.d-1
!      radius(3) = 5.d-2
!      phi(1) = 0.d0
!      phi(2) = 0.d0
!      phi(3) = pi/4.d0

      radius(1) = 2.d-1
      radius(2) = 1.d-1
      centerx(1) = 0.5d0
      centerx(2) = 0.5d0 + 0.3d0 - 0.006d0
      centery(1) = 0.1d0
      centery(2) = 0.1d0

c      radius(3) = 5.d-2
      phi(1) = 0.d0
      phi(2) = 0.d0
c      phi(3) = pi/4.d0
      
c      Nphi = 128
c      dphi = twopi/Nphi
c      open(unit=21,file='output/dragx.dat')
c      open(unit=22,file='output/dragy.dat')

c      do kk = 1,Nphi
c        phi(1) = dble(kk-1)*dphi

      do j = 1,nbodies
        do k = 1,ninner
          theta = dble(k-1)*dtheta
          var_rad = radius(j)
          x((j-1)*ninner+k) = centerx(j) + var_rad*
     $          (dcos(phi(j))*dcos(theta) + 
     $           2.d0*dsin(phi(j))*dsin(theta))
          y((j-1)*ninner+k) = centery(j) + var_rad*
     $          (-dsin(phi(j))*dcos(theta) + 
     $           2.d0*dcos(phi(j))*dsin(theta))
        enddo
      enddo


      xmin = -8.d-1
      xmax = 8.d-1
      ymin = -0.99d0
      ymax = 0.99d0
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

c      icount = 0
c      do j = 1,30
c        var_rad = radius(1)*(1.04d0 + 1.d-2*j)
c        do k = 1,ninner
c          icount = icount + 1
c          theta = dble(k-1)*dtheta
c          xtar(icount) = centerx(1) + var_rad*
c     $        (dcos(phi(1))*dcos(theta) + 
c     $        2.d0*dsin(phi(1))*dsin(theta))
c          ytar(icount) = centery(1) + var_rad*
c     $        (-dsin(phi(1))*dcos(theta) + 
c     $        2.d0*dcos(phi(1))*dsin(theta))
c        enddo
c      enddo

      do j = 1,2*ninner*nbodies + 3*nbodies + 2*nouter
        den(j) = 0.d0
      enddo

      ifmm = 2
      call stokesSolver(ninner,nbodies,nouter,ifmm,x,y,den,iter)
c     pass in number of points and x and y coordinates and return the
c     density function on the boundary
      print *, 'Solver'
c      call classifyPoints(ninner,nbodies,x,y,
c     $      ntargets*ntargets,xtar,ytar,iside,inear)
c
c      call computeVorticityTargets(ninner,nbodies,nouter,
c     $      x,y,den,ntargets*ntargets,xtar,ytar,vort_tar)
c
c      call computeVelocityPressureTargets(ninner,nbodies,nouter,
c     $    x,y,den,ntargets*ntargets,xtar,ytar,utar,vtar,press_tar)
c
      call computeQoiTargetsBaryCen(ninner,nbodies,nouter,x,y,den,
     $  ntargets*ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
      print *, 'computeQoiTargets'
      
      call computeShearStress(ninner,nbodies,nouter,x,y,den,
     $    shear_stress)
c     pass in the density function and return the shear_stress
      print *, 'computeShearStress'
      call computePressure(ninner,nbodies,nouter,x,y,den,pressure)
c     pass in the density function and return the pressure
      print *, 'computePressure'
      call computeDrag(ninner,nbodies,x,y,
     $      shear_stress,pressure,drag)
c     pass in the shear_stress and pressure and return the drag
c      write(21,1000) drag(1)
c      write(22,1000) drag(2)
c      enddo
c      close(unit=21)
c      close(unit=22)
      print *, 'computeDrag' 
      
      
      if (1 .eq. 1) then
        do j = 1,nbodies        
         do k = 1,ninner
          theta = dble(k-1)*dtheta
          var_rad = radius(j) + 0.000001d0
          xtar_test((j-1)*ninner+k) = centerx(j) + var_rad*
     $          (dcos(phi(j))*dcos(theta) + 
     $           2.d0*dsin(phi(j))*dsin(theta))
          ytar_test((j-1)*ninner+k) = centery(j) + var_rad*
     $          (-dsin(phi(j))*dcos(theta) + 
     $           2.d0*dcos(phi(j))*dsin(theta))
         enddo
       enddo
      call computeQoiTargetsBaryCen(ninner,nbodies,nouter,x,y,den,
     $  ninner*nbodies,xtar_test,ytar_test,utar_test,vtar_test,
     $  press_tar_test,vort_tar_test) 
     
       open(unit=10,file='output/shear_stress_test.dat')     
       do k = 1,ninner*nbodies
         write(10,1000) vort_tar_test(k)
       enddo      
       close(unit=10)
      endif

      open(unit=1,file='output/den.dat')
      open(unit=2,file='output/shear_stress.dat')
      open(unit=3,file='output/pressure.dat')
      open(unit=4,file='output/drag.dat')
      open(unit=8,file='output/x.dat')
      open(unit=9,file='output/y.dat')
 
c      open(unit=1,file='output/trap_den.dat') 
c      open(unit=2,file='output/trap_shear_stress.dat') 
c      open(unit=3,file='output/trap_pressure.dat') 
c      open(unit=4,file='output/trap_drag.dat')      
c      open(unit=8,file='output/trap_x.dat')      
c      open(unit=9,file='output/trap_y.dat')      
      
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
      open(unit=9,file='output/vort_tar.dat')
      
c      open(unit=1,file='output/trap_xtar.dat')      
c      open(unit=2,file='output/trap_ytar.dat')      
c      open(unit=3,file='output/trap_utar.dat')      
c      open(unit=4,file='output/trap_vtar.dat')      
c      open(unit=8,file='output/trap_press_tar.dat')      
c      open(unit=9,file='output/trap_vort_tar.dat')      

      do k = 1,ntargets**2
        write(1,1000) xtar(k)
        write(2,1000) ytar(k)
        write(3,1000) utar(k)
        write(4,1000) vtar(k)
        write(8,1000) press_tar(k)
        write(9,1000) vort_tar(k)
      enddo
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=8)
      close(unit=9)

c      write(6,1010) phi(1),drag(1),drag(2)

 1000 format(E25.16)
 1010 format(ES20.4)



      end

