      program stokesDriver
      implicit real*8 (a-h,o-z)

c      parameter (ninner = 512)
      parameter (ninner = 2**13)
      parameter (nbodies = 1)
      parameter (maxtargets = 20000)

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      dimension centerx(20),centery(20)
      dimension radius(20),phi(20)
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension shear_stress(ninner*nbodies)
      dimension xtar(maxtargets),ytar(maxtargets)
      dimension utar(maxtargets),vtar(maxtargets)
      dimension press_tar(maxtargets)
      dimension x0(maxtargets),y0(maxtargets)


      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(ninner)

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)

      call srand(seed(1))
c      print*,rand(),rand(),rand(),rand()

c      centerx(1) = -6.d-1
c      centerx(2) = -3.d-1
c      centerx(3) = -2.d-1
c      centerx(4) = 0.d-1
c      centerx(5) = 1.d-1
c      centerx(6) = 2.d-1
c      centerx(7) = 4.d-1
c      centerx(8) = 5.d-1
c      centerx(9) = 1.d-1
c      centerx(10) = -2.d-1
c      centerx(11) = 8.d-1
c      centerx(12) = -4.d-1
c      centerx(13) = -5.d-1
c      centerx(14) = 4.d-1
c      centerx(15) = 1.5d-1
c      centerx(16) = 5.d-1
c      centerx(17) = 7.d-1
c      centerx(18) = -6.d-1
c      centerx(19) = 0.d-1
c      centerx(20) = 4.d-1
c
c      centery(1) = 2.d-1
c      centery(2) = -1.d-1
c      centery(3) = 3.d-1
c      centery(4) = -2.d-1
c      centery(5) = 6.d-1
c      centery(6) = -4.d-1
c      centery(7) = 3.d-1
c      centery(8) = -4.d-1
c      centery(9) = 1.d-1
c      centery(10) = -4.d-1
c      centery(11) = 2.d-1
c      centery(12) = 6.d-1
c      centery(13) = -5.d-1
c      centery(14) = -1.d-1
c      centery(15) = 3.5d-1
c      centery(16) = 6.d-1
c      centery(17) = -1.d-1
c      centery(18) = -1.d-1
c      centery(19) = -7.d-1
c      centery(20) = -6.5d-1
c
c      do k = 1,20
c        radius(k) = 8.d-2 + 2.d-2*(rand()-5.d-1)
c        phi(k) = 0.d0*twopi*rand()
c        centerx(k) = centerx(k) + 2.d-2*(rand()-5.d-1)
c        centery(k) = centery(k) + 2.d-2*(rand()-5.d-1)
c      enddo
c
cc      centerx(1) = -0.5d0
cc      centery(1) = 0.0d0
cc      centerx(2) = 0.5d0
cc      centery(2) = -0.0d0
ccc      centerx(3) = -0.0d0
ccc      centery(3) = -0.2d0
cc      radius(1) = 2.d-1
cc      radius(2) = 2.d-1
ccc      radius(3) = 1.d-1
cc      phi(1) = 0.d0
cc      phi(2) = 0.d0
ccc      phi(3) = 0.d0

      centerx(1) = 0.d0
      centery(1) = 0.d0
      phi(1) = 0.d0
c      radius(1) = 8.d-1

      nx = 100 
      ny = 41
      ntargets = nx*ny
      xmin = -9.0d0
      xmax = +9.0d0
      ymin = -8.0d-1
      ymax = +8.0d-1
      dx = (xmax - xmin)/(nx-1)
      dy = (ymax - ymin)/(ny-1)
      itar = 0
      do j = 1,nx
        do k = 1,ny
          itar = itar + 1
          xtar(itar) = xmin + dble(j-1)*dx
          ytar(itar) = ymin + dble(k-1)*dy
        enddo
      enddo
      
c      ntargets = 80
c      nhalf = ntargets/2
c      ymax = 8.0d-1
c      dy = 2*ymax/(nhalf - 1)
c      do k = 1,nhalf
c        xtar(k) = -9.7d0
c        ytar(k) = -1.d0*ymax + dble(k-1)*dy
c      enddo
c      do k = nhalf+1,2*nhalf
c        xtar(k) = 9.7d0
c        ytar(k) = -1.d0*ymax + dble(k-nhalf-1)*dy
c      enddo

      open(unit=11,file='output/radii.dat')
      open(unit=12,file='output/press_drop.dat')

      smoothOrder = 8.d0
      dr = 1.d-2
c      radius(1) = 0.d0
      radius(1) = 5.d-1 - dr
      do i=1,1
        radius(1) = radius(1) + dr

        do j = 1,nbodies
          do k = 1,ninner
            theta = dble(k-1)*dtheta
c            var_rad = radius(j)*(1.d0 + 2.d-1*dcos(5*theta))
            var_rad = radius(j)
            var_rad = var_rad*(dcos(theta)**smoothOrder + 
     $          dsin(theta)**smoothOrder)**(-1.d0/smoothOrder)
            x((j-1)*ninner+k) = centerx(j) + var_rad*
     $          (dcos(phi(j))*dcos(theta) + dsin(phi(j))*dsin(theta))
            y((j-1)*ninner+k) = centery(j) + var_rad*
     $          (-dsin(phi(j))*dcos(theta) + dcos(phi(j))*dsin(theta))
c            x((j-1)*ninner+k) = var_rad*dcos(theta) + centerx(j)
c            y((j-1)*ninner+k) = var_rad*dsin(theta) + centery(j) 
          enddo
        enddo


        ifmm = 1
        call stokesSolver(ninner,nbodies,ntargets,ifmm,x,y,
     $        xtar,ytar,utar,vtar,press_tar,shear_stress)
c       pass in number of points and x and y coordinates and return the
c       shear stress on the boundary
        ave_left = 0.d0
        ave_right = 0.d0
        do k=1,nhalf
          ave_left =  ave_left + press_tar(k)
        enddo
        ave_left = ave_left/dble(nhalf)
        do k=nhalf+1,2*nhalf
          ave_right = ave_right + press_tar(k)
        enddo
        ave_right = ave_right/dble(nhalf)
        press_drop = ave_left - ave_right

        write(11,1000) radius(1)
        write(12,1000) press_drop/(xtar(nhalf+1)-xtar(1))
        print*,radius(1)


c        ntime = 2
c        dt = 1.d-1
c        x0(1) = -8.d-1
c        y0(1) = 5.d-1
cc       initial condition
c        call stokesTimeStepper(ninner,nbodies,x,y,ntime,dt,
c     $      x0,y0)



      enddo

      close(unit=11)
      close(unit=12)
 1000 format(E25.16)

      end

