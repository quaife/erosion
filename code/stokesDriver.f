      program stokesDriver
      implicit real*8 (a-h,o-z)

      parameter (ninner = 512)
      parameter (nbodies = 20)

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      dimension centerx(20),centery(20)
      dimension radius(20),phi(20)
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension shear_stress(ninner*nbodies)

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)

      call srand(seed(1))
c      print*,rand(),rand(),rand(),rand()

      centerx(1) = -6.d-1
      centerx(2) = -3.d-1
      centerx(3) = -2.d-1
      centerx(4) = 0.d-1
      centerx(5) = 1.d-1
      centerx(6) = 2.d-1
      centerx(7) = 4.d-1
      centerx(8) = 5.d-1
      centerx(9) = 1.d-1
      centerx(10) = -2.d-1
      centerx(11) = 8.d-1
      centerx(12) = -4.d-1
      centerx(13) = -5.d-1
      centerx(14) = 4.d-1
      centerx(15) = 1.5d-1
      centerx(16) = 5.d-1
      centerx(17) = 7.d-1
      centerx(18) = -6.d-1
      centerx(19) = 0.d-1
      centerx(20) = 4.d-1

      centery(1) = 2.d-1
      centery(2) = -1.d-1
      centery(3) = 3.d-1
      centery(4) = -2.d-1
      centery(5) = 6.d-1
      centery(6) = -4.d-1
      centery(7) = 2.d-1
      centery(8) = -4.d-1
      centery(9) = 1.d-1
      centery(10) = -4.d-1
      centery(11) = 2.d-1
      centery(12) = 6.d-1
      centery(13) = -5.d-1
      centery(14) = -1.d-1
      centery(15) = 3.5d-1
      centery(16) = 6.d-1
      centery(17) = -1.d-1
      centery(18) = -1.d-1
      centery(19) = -7.d-1
      centery(20) = -6.5d-1

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(ninner)
      do k = 1,20
        radius(k) = 8.d-2 + 2.d-2*(rand()-5.d-1)
        phi(k) = twopi*rand()
        centerx(k) = centerx(k) + 2.d-2*(rand()-5.d-1)
        centery(k) = centery(k) + 2.d-2*(rand()-5.d-1)
      enddo

      do j = 1,nbodies
        do k = 1,ninner
          theta = dble(k-1)*dtheta
          var_rad = radius(j)*(1.d0 + 2.d-1*dcos(5*theta))
c          var_rad = radius(j)
          x((j-1)*ninner+k) = centerx(j) + var_rad*
     $        (dcos(phi(j))*dcos(theta) + dsin(phi(j))*dsin(theta))
          y((j-1)*ninner+k) = centery(j) + var_rad*
     $        (-dsin(phi(j))*dcos(theta) + dcos(phi(j))*dsin(theta))
c          x((j-1)*ninner+k) = var_rad*dcos(theta) + centerx(j)
c          y((j-1)*ninner+k) = var_rad*dsin(theta) + centery(j) 
        enddo
      enddo
      

c      twopi = 8.d0*datan(1.d0)
c      dtheta = twopi/dble(ninner)
c      do j = 1,nbodies
c        if (j .eq. 1) then
c          rad = 2.d-1
c        elseif (j .eq. 2) then
c          rad = 1.5d-1
c        else
c          rad = 1.d-1
c        endif
c        do k = 1,ninner
c          theta = dble(k-1)*dtheta
cc          var_rad = rad
c          var_rad = rad*(1.d0 + 2.d-1*dcos(5*theta))
c          x((j-1)*ninner+k) = var_rad*dcos(theta) + 6.d-1*(j-2)
c          y((j-1)*ninner+k) = var_rad*dsin(theta) + 3.d-1*(j-2)
c        enddo
c      enddo
cc     parameterize from left most point and proceed up the top of the
cc     curve, throught the back point, and back to the leading
cc     singularity point

      call stokesSolver(ninner,nbodies,x,y,shear_stress)
c     pass in number of points and x and y coordinates and return the
c     shear stress on the boundary

      end

