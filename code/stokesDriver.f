      program stokesDriver
      implicit real*8 (a-h,o-z)

      parameter (ninner = 128)
      parameter (nbodies = 20)
      parameter (ntargets = 100)

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      dimension centerx(20),centery(20)
      dimension radius(20),phi(20)
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension shear_stress(ninner*nbodies)
      dimension xtar(ntargets),ytar(ntargets)
      dimension utar(ntargets),vtar(ntargets)
      dimension press_tar(ntargets)

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
      centery(7) = 3.d-1
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
        phi(k) = 0.d0*twopi*rand()
        centerx(k) = centerx(k) + 2.d-2*(rand()-5.d-1)
        centery(k) = centery(k) + 2.d-2*(rand()-5.d-1)
      enddo

c      centerx(1) = -0.5d0
c      centery(1) = 0.0d0
c      centerx(2) = 0.5d0
c      centery(2) = -0.0d0
cc      centerx(3) = -0.0d0
cc      centery(3) = -0.2d0
c      radius(1) = 2.d-1
c      radius(2) = 2.d-1
cc      radius(3) = 1.d-1
c      phi(1) = 0.d0
c      phi(2) = 0.d0
cc      phi(3) = 0.d0

      do j = 1,nbodies
        do k = 1,ninner
          theta = dble(k-1)*dtheta
c          var_rad = radius(j)*(1.d0 + 2.d-1*dcos(5*theta))
          var_rad = radius(j)
          x((j-1)*ninner+k) = centerx(j) + var_rad*
     $        (dcos(phi(j))*dcos(theta) + dsin(phi(j))*dsin(theta))
          y((j-1)*ninner+k) = centery(j) + var_rad*
     $        (-dsin(phi(j))*dcos(theta) + dcos(phi(j))*dsin(theta))
c          x((j-1)*ninner+k) = var_rad*dcos(theta) + centerx(j)
c          y((j-1)*ninner+k) = var_rad*dsin(theta) + centery(j) 
        enddo
      enddo
c
c
c      nx = 20 
c      ny = 50
c
c      xmin = -2.7d0
c      xmax = +2.7d0
c      ymin = -8.0d-1
c      ymax = +8.0d-1
c      dx = (xmax - xmin)/(nx-1)
c      dy = (ymax - ymin)/(ny-1)
c      itar = 0
c      do j = 1,nx
c        do k = 1,ny
c          itar = itar + 1
c          xtar(itar) = xmin + dble(j-1)*dx
c          ytar(itar) = ymin + dble(k-1)*dy
c        enddo
c      enddo

      
      nhalf = ntargets/2
      ymax = 8.0d-1
      dy = 2*ymax/(nhalf - 1)
      do k = 1,nhalf
        xtar(k) = -2.7d0
        ytar(k) = -1.d0*ymax + dble(k-1)*dy
      enddo
      do k = nhalf+1,2*nhalf
        xtar(k) = 2.7d0
        ytar(k) = -1.d0*ymax + dble(k-nhalf-1)*dy
      enddo

      call stokesSolver(ninner,nbodies,ntargets,x,y,
     $      xtar,ytar,utar,vtar,press_tar,shear_stress)
c     pass in number of points and x and y coordinates and return the
c     shear stress on the boundary

      end

