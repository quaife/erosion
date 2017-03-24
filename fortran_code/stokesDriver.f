      program stokesDriver
      implicit real*8 (a-h,o-z)

c      parameter (ninner = 512)
      parameter (ninner = 32)
      parameter (nbodies = 1)
      parameter (nouter = 2**12)

      dimension centerx(20),centery(20)
      dimension radius(20),phi(20)
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension den(2*ninner*nbodies + 3*nbodies + 2*nouter)
      dimension shear_stress(ninner*nbodies)
      dimension pressure(ninner*nbodies)
      dimension drag(2*nbodies)

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(ninner)

      centerx(1) = -0.5d0
      centery(1) = 0.0d0
      centerx(2) = 0.5d0
      centery(2) = -0.0d0
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
          y((j-1)*ninner+k) = centery(j) + var_rad*
     $          (-dsin(phi(j))*dcos(theta) + dcos(phi(j))*dsin(theta))
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


      open(unit=1,file='output/den.dat')
      open(unit=2,file='output/shear_stress.dat')
      open(unit=3,file='output/pressure.dat')
      open(unit=4,file='output/drag.dat')
      do k = 1,2*ninner*nbodies + 2*nouter + 3*nbodies
        write(1,1000) den(k)
      enddo
      do k = 1,ninner*nbodies
        write(2,1000) shear_stress(k)
        write(3,1000) pressure(k)
      enddo
      do k = 1,2*nbodies
        write(4,1000) drag(k)
      enddo
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)



 1000 format(E25.16)



      end

