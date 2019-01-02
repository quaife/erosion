      subroutine stokesSolver(nninner,nnbodies,nnouter,
     $     ifmm,xx,yy,den,iter)
c     Input x and y coordinates and return the density function on the
c     boundary.  Outer wall is used to enclose the inner boundary so
c     that Stokes paradox is avoided
      implicit real*8 (a-h,o-z)

      dimension xx(nninner*nnbodies),yy(nninner*nnbodies)
c     x and y coordinates of obstacle

      parameter (nmax = 2**17)
      parameter (maxbodies = 50)
      parameter (ntargets = 2500)
c     max points on the boundary of the obstacle      
      parameter (maxl = 2000, liwork = 30)
      parameter (lrwork = 10 + nmax*(maxl+6) + 
     $     maxl*(maxl+3))
c     maximum size of workspaces for GMRES
c     maxl is the maximum number of GMRES steps

      dimension x(nmax),y(nmax)
      dimension centerx(maxbodies),centery(maxbodies)
      dimension px(nmax),py(nmax)
c     x and y coordinates of the normal of the obstacle
      dimension cur(nmax),speed(nmax)
c     Jacobian and curvature of the geometry

c      dimension xtar(ntargets),ytar(ntargets)
cc     Target locations
c      dimension utar(ntargets),vtar(ntargets)
c      dimension press_tar(ntargets)
cc     Velocity and pressure at target locations

      dimension xouter(nmax),youter(nmax)
c     x and y coordinates of confining wall
      dimension px0(nmax), py0(nmax)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nmax), speed0(nmax)
c     Jacobian and curvature of the confining wall

      dimension gmwork(lrwork), igwork(liwork)
c     workspaces for GMRES

      dimension rhs(nmax)
c     right hand side
      dimension den(nmax)

cc     x and y coordinates of the density function
c      dimension E11(nninner*nnbodies)
c      dimension E12(nninner*nnbodies)
c      dimension E22(nninner*nnbodies)
cc     components of the deformation gradient on the obstacle
c      dimension shear_stress(nninner*nnbodies)
cc     shear stress on the obstacle

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter
c     global variables for the geometry which are needed by the external
c     matvec routine in gmres

      ninner = nninner
      nbodies = nnbodies
      nouter = nnouter
c     can't declare input variables into a common field, so need new
c     variable name for x,y,ninner,nbodies
      call inner_geometry(ninner,nbodies,x,y,xx,yy,px,py,cur,speed,
     $    centerx,centery)
c
      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     load geometry of initial shape


c      open(unit=1,file='output/xinner.dat')
c      open(unit=2,file='output/yinner.dat')
c      open(unit=3,file='output/xouter.dat')
c      open(unit=4,file='output/youter.dat')
c      do k = 1,ninner*nbodies
c        write(1,1000) x(k)
c        write(2,1000) y(k)
c      enddo
c      do k = 1,nouter
c        write(3,1000) xouter(k)
c        write(4,1000) youter(k)
c      enddo
c      close(unit=1)
c      close(unit=2)
c      close(unit=3)
c      close(unit=4)

      call bd_condition(ninner,nbodies,x,y,nouter,xouter,youter,rhs)
c     load boundary condition

      call solveBIE(ninner,nbodies,nouter,den,rhs,
     $      gmwork,lrwork,igwork,liwork,maxl,ifmm,ibary,iter)
c     solve for the density function with GMRES

c      call eval_velocity_targets(ninner,nbodies,
c     $    x,y,centerx,centery,px,py,speed,
c     $    nouter,xouter,youter,px0,py0,speed0,den,
c     $    ntargets,xtar,ytar,utar,vtar,press_tar)
c
c      open(unit=1,file='output/xtar.dat')
c      open(unit=2,file='output/ytar.dat')
c      open(unit=3,file='output/utar.dat')
c      open(unit=4,file='output/vtar.dat')
c      open(unit=5,file='output/presstar.dat')
c
c      do k = 1,ntargets
c        write(1,1000) xtar(k)
c        write(2,1000) ytar(k)
c        write(3,1000) utar(k)
c        write(4,1000) vtar(k)
c        write(5,1000) press_tar(k)
c      enddo
c
c      close(unit=1)
c      close(unit=2)
c      close(unit=3)
c      close(unit=4)
c      close(unit=5)
      

c 1000 format(E25.16)
      end



c***********************************************************************
      subroutine inner_geometry(ninner,nbodies,x,y,
     $    xx,yy,px,py,cur,speed,centerx,centery)
      implicit real*8 (a-h,o-z)

      dimension xx(ninner*nbodies),yy(ninner*nbodies)
      dimension x(*),y(*)
      dimension px(*), py(*)
      dimension cur(*), speed(*)
      dimension centerx(*),centery(*)

      do k = 1,ninner*nbodies
        x(k) = xx(k)
        y(k) = yy(k)
      enddo

      do ibod = 1,nbodies
        centerx(ibod) = 0.d0
        centery(ibod) = 0.d0
        do k = 1,ninner
          centerx(ibod) = centerx(ibod) + x((ibod-1)*ninner+k)
          centery(ibod) = centery(ibod) + y((ibod-1)*ninner+k)
        enddo
        centerx(ibod) = centerx(ibod)/dble(ninner)
        centery(ibod) = centery(ibod)/dble(ninner)
      enddo
c     Find the center of each body for stokeslets and rotlets

      call compute_derivs(ninner,nbodies,x,y,px,py,cur,speed)
c     spectrally compute normal, curvature, and Jacobian of curve

      end

c***********************************************************************
      subroutine outer_geometry(nouter,x,y,px,py,cur,speed)
c     Load the outer shape of the geometry
      implicit real*8 (a-h,o-z)

      dimension x(*),y(*)
      dimension px(*), py(*)
      dimension cur(*), speed(*)

      complex *16 eye

      eye = (0.d0,1.d0)
      twopi = 8.d0*datan(1.d0)

      dtheta = twopi/dble(nouter)

      do k = 1,nouter
        theta = dble(k-1)*dtheta
        smoothOrder = 8.d0
        r = (dcos(theta)**smoothOrder+ dsin(theta)**smoothOrder)**
     $       (-1.d0/smoothOrder)
        x(k) = 3.d0*r*dcos(theta)
        y(k) = r*dsin(theta)
      enddo
c     outer geometry

      call compute_derivs(nouter,1,x,y,px,py,cur,speed)
      do k = 1,nouter
        px(k) = -1.d0*px(k)
        py(k) = -1.d0*py(k)
        cur(k) = -1.d0*cur(k)
      enddo
c     normal and curvatures need to point outwards



      end
        

c***********************************************************************
      subroutine bd_condition(ninner,nbodies,x,y,
     $    nouter,xouter,youter,rhs)
c     Boundary condition coming from the far field flow
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension xouter(nouter),youter(nouter)
      dimension den(2*nouter + 2*ninner*nbodies+3)
      dimension rhs(2*nouter + 2*ninner*nbodies+3)

      twopi = 8.d0*datan(1.d0)

      if (0 .eq. 1) then
        sto1 = 1.d0
        sto2 = 2.d0
        rot = 3.d0

        cx = 0.0d0
        cy = 0.0d0
        do k = 1,nouter
          rx = xouter(k) - cx
          ry = youter(k) - cy
          rdots = rx*sto1 + ry*sto2
          rho2 = rx**2.d0 + ry**2.d0

          rhs(k) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx) +
     $      rot*ry/rho2
          rhs(k+nouter) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry) -
     $      rot*rx/rho2
        enddo

        do j = 1,nbodies
          do k=1,ninner
            rx = x((j-1)*ninner+k) - cx
            ry = y((j-1)*ninner+k) - cy
            
            rdots = rx*sto1 + ry*sto2
            rho2 = rx**2.d0 + ry**2.d0

            rhs(2*nouter+(j-1)*2*ninner+k) = 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx) +
     $        rot*ry/rho2
            rhs(2*nouter+(j-1)*2*ninner+k+ninner) = 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry) -
     $        rot*rx/rho2
          enddo
        enddo
      endif
c     velocity coming from the r \otimes r part of the Stokes flow

      if (0 .eq. 1) then
        do k = 1,nouter
          call bgFlow(1,xouter(k),youter(k),rhs(k),rhs(k+nouter))
        enddo

        do j = 1,nbodies
          do k = 1,ninner
            call bgFlow(1,x((j-1)*ninner+k),y((j-1)*ninner+k),
     $          rhs(2*nouter+(j-1)*2*ninner+k),
     $          rhs(2*nouter+(j-1)*2*ninner+k+ninner))
          enddo
        enddo
      endif
c     generic boundary condition

      if (1 .eq. 1) then
        do k = 1,nouter
          call bgFlow(1,xouter(k),youter(k),rhs(k),rhs(k+nouter))
        enddo
        do j = 1,nbodies
          do k = 1,ninner
            rhs(2*nouter + (j-1)*2*ninner + k) = 0.d0
            rhs(2*nouter + (j-1)*2*ninner + k + ninner) = 0.d0
          enddo
        enddo
c        do k = 1,ninner*nbodies
c          rhs(k+2*nouter) = 0.d0
c          rhs(k+2*nouter+ninner) = 0.d0
c        enddo
      endif
c     constant background flow on outer walls on no slip on obstacle

      if (0 .eq. 1) then
        do k = 1,nouter
          rho2 = (xouter(k)-2.d0)**2.d0 + (youter(k)+2.d0)**2.d0
          rhs(k) = ((youter(k)+2.d0)**2.d0 - (xouter(k)-2.d0)**2.d0)/
     $        rho2**2.d0
          rhs(k + nouter) = -2.d0*(xouter(k)-2.d0)*(youter(k)+2.d0)/
     $        rho2**2.d0
        enddo

        do k = 1,ninner*nbodies
          rho2 = (x(k)-2.d0)**2.d0 + (y(k)+2.d0)**2.d0
          rhs(2*nouter + k) = ((y(k)+2.d0)**2.d0 - (x(k)-2.d0)**2.d0)/
     $        rho2**2.d0
          rhs(2*nouter + k + ninner) = -2.d0*(x(k)-2.d0)*(y(k)+2.d0)/
     $        rho2**2.d0
        enddo
      endif


      do k = 1,nbodies
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+1) = 0.d0
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+2) = 0.d0
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+3) = 0.d0
      enddo
c     extra terms need for rotlet and stokeslet conditions

      end


c***********************************************************************
      subroutine fourierUpsample(n,nrate,zIn,zOut)
      implicit real*8 (a-h,o-z)

      complex *16 zIn(n)
      complex *16 zOut(n*nrate)

      dimension wsave(4*n*nrate + 15)

      call DCFFTI(n,wsave)
      call DCFFTF(n,zIn,wsave)
c     Take forward FFT

      do k = 1,n*nrate
        zOut(k) = 0.d0
      enddo

      do k = 1,n/2
        zOut(k) = zIn(k)/dble(n)
        zOut(n*nrate - k + 1) = zIn(n - k + 1)/dble(n)
      enddo

      call DCFFTI(n*nrate,wsave)
      call DCFFTB(n*nrate,zOut,wsave)
c     Take backward FFT



      end

c***********************************************************************
      subroutine fourierDiff(n,den,wsave)
      implicit real*8 (a-h,o-z)

      complex *16 den(n)
      real *8 wsave(4*n+15)
      complex *16 eye

      dimension modes(n)
      real *8 work(2*n)
      lenwork = 2*n

      eye = (0.d0,1.d0)

      do k = 1,n/2
        modes(k) = k-1
      enddo
      do k = n/2+1,n
        modes(k) = -n+k-1
      enddo

      call DCFFTI(n,wsave)
c     Initialize variable for FFT

      call DCFFTF(n,den,wsave)
c     Take forward FFT

      do k = 1,n/2
        den(k) = den(k)*eye*dble(k-1)/dble(n)
      enddo
c     Multiply by 1i*modes for positive frequencies
      den(n/2+1) = 0.d0
c     kill the Nyquist frequency
      do k=n/2+2,n
        den(k) = den(k)*eye*dble(-n+k-1)/dble(n)
      enddo
c     Multiply by 1i*modes for negative frequencies

      call DCFFTB(n,den,wsave)
c     Take backward FFT

      end


c***********************************************************************
      subroutine bgFlow(ninner,x,y,u,v)
      implicit real*8 (a-h,o-z)

      dimension x(ninner),y(ninner)
      dimension u(ninner),v(ninner)

      do k=1,ninner
c        u(k) = 1.d0
c        v(k) = 0.d0
cc       constant flow

c        u(k) = x(k)
c        v(k) = -y(k)
cc       extesional flow

        scal = 1.0d0
        u(k) = scal*(1.0d0 + y(k))*(1.0d0 - y(k))
        v(k) = 0.d0
c       pipe flow

c        u(k) = y(k)
c        v(k) = 0.d0
cc       shear flow
    
c        u(k) = y(k)/(x(k)**2.d0 + y(k)**2.d0)
c        v(k) = -x(k)/(x(k)**2.d0 + y(k)**2.d0)
cc       single Rotlet
      enddo

      end 


c***********************************************************************
      subroutine compute_derivs(n,nbodies,x,y,px,py,cur,speed)
      implicit real*8 (a-h,o-z)

      dimension x(n*nbodies),y(n*nbodies)
      dimension px(n*nbodies),py(n*nbodies)
      dimension cur(n*nbodies),speed(n*nbodies)

      real *8 Dx(n),Dy(n)
      real *8 DDx(n),DDy(n)
      complex *16 zden1(n),zden2(n)
      real*8 wsave(4*n+15)
      complex *16 eye

      eye = (0.d0,1.d0)

      do j = 1,nbodies
        do k = 1,n
          zden1(k) = x((j-1)*n+k) + eye*y((j-1)*n+k)
        enddo
        call fourierDiff(n,zden1,wsave)
        zden2 = zden1
        call fourierDiff(n,zden2,wsave)
c       real part of zden is the first derivative of the x component of
c       the position, imaginary part of zden is the first derivative of
c       the y component of the position

        Dx = dreal(zden1)
        Dy = dimag(zden1)
        DDx = dreal(zden2)
        DDy = dimag(zden2)
        do k=1,n
          speed((j-1)*n+k) = dsqrt(Dx(k)**2.d0 + Dy(k)**2.d0)
          px((j-1)*n+k) = -Dy(k)/speed((j-1)*n+k)
          py((j-1)*n+k) = Dx(k)/speed((j-1)*n+k)
          cur((j-1)*n+k) = (Dy(k)*DDx(k) - Dx(k)*DDy(k))/
     $        speed((j-1)*n+k)**3.d0
        enddo
      enddo
        

      end


c***********************************************************************
      subroutine solveBIE(ninner,nbodies,nouter,den,rhs,
     $    gmwork,lrwork,igwork,liwork,maxl,ifmm,ibary,iter)
c     Solve the boundary integral equation with GMRES
      implicit real*8 (a-h,o-z)

      external msolve_DLP,matvec_DLP,matvec_DLP_fmm

      dimension den(2*nouter+2*ninner*nbodies+3*nbodies)
      dimension rhs(2*nouter+2*ninner*nbodies+3*nbodies)
c     leave room for stokeslets and rotlets at the end of den
c      dimension vel1(2*nouter+2*ninner*nbodies+3*nbodies)
c      dimension vel2(2*nouter+2*ninner*nbodies+3*nbodies)

      dimension gmwork(lrwork),igwork(liwork)
c     gmres workspaces

      dimension iwork(3)
c     workspace integer array that'll be passed to preconditioner msolve

      itol = 0
      tol = 1.d-8
      isym = 0
      igwork(1) = maxl
      do i=2,liwork
        igwork(i) = 0
      enddo
c     paramters for DGMRES

      igwork(4) = 0
c     no preconditioner
c      igwork(4) = -1
cc     left precondition
c      igwork(4) = 1
cc     right precondition

      igwork(5) = -1
c     restart flag

      iwork(1) = ninner
      iwork(2) = nbodies
      iwork(3) = nouter

cc     DEBUG
c      twopi = 8.d0*datan(1.d0)
c      dtheta = twopi/dble(ninner)
c      do k = 1,ninner
c        theta = dble(k-1)*dtheta
c        rhs(2*nouter + k) = dexp(dcos(theta))
c        rhs(2*nouter + k + ninner) = dexp(dsin(theta))
c      enddo

c      do k = 1,2*nouter+2*ninner*nbodies+3*nbodies
c        den(k) = rhs(k)
c      enddo
c     initial guess

cc      DEBUG
c      ntotal = 2*nouter + 2*ninner*nbodies + 3*nbodies 
c      call matvec_DLP(ntotal,den,vel1,nelt,ia,ja,a,isym)
c      call matvec_DLP_fmm(ntotal,den,vel2,nelt,ia,ja,a,isym)
c      do k = 1,ntotal
cc        if (abs(vel1(k) - vel2(k)) .ge. 1.d-10) then
c          print*,k,vel1(k)-vel2(k)
cc        endif
c      enddo

      if (ibary .eq. 1) then
        print*,'USING Barycentric'
        call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
     $      nelt,ia,ja,a,isym,
     $      matvec_DLP_bary,msolve_DLP,itol,tol,itmax,iter,err,
     $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c       Use Barycentric
      else

        if (ifmm .eq. 1) then
          print*,'USING FMM'
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
       $      nelt,ia,ja,a,isym,
       $      matvec_DLP_fmm,msolve_DLP,itol,tol,itmax,iter,err,
       $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c         USE FMM
        else
          print*,'USING DIRECT'
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
       $      nelt,ia,ja,a,isym,
       $      matvec_DLP,msolve_DLP,itol,tol,itmax,iter,err,
       $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c          DON'T USE FMM
        endif
c       use GMRES to find the density function of the double layer
c       potential
      endif

      end


c***********************************************************************
      subroutine matvec_DLP(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 50)

      dimension den(ntotal)
      dimension vel(ntotal)

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter

      dimension x(nmax),y(nmax)
      dimension px(nmax),py(nmax)
      dimension cur(nmax),speed(nmax)
      dimension centerx(maxbodies),centery(maxbodies)

      dimension xouter(nmax),youter(nmax)
      dimension px0(nmax),py0(nmax)
      dimension cur0(nmax),speed0(nmax)

      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))
      dimension ux(max(ninner,nouter)),uy(max(ninner,nouter))

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do k = 1,2*nouter + 2*ninner*nbodies + 3*nbodies
        vel(k) = 0.d0
      enddo
c     initialize velocity to 0

c     START OF SOURCE POINTS == OUTER WALL
      do k = 1,nouter
        denx(k) = den(k)
        deny(k) = den(k + nouter)
      enddo
c     density function due to outer wall

c     START OF TARGET POINTS == OUTER WALL
c     loop over target points
      do k=1,nouter
        ux(k) = 0.d0
        uy(k) = 0.d0

c       loop over source points
        do j=1,nouter
          if (j .ne. k) then
            rx = xouter(k) - xouter(j)
            ry = youter(k) - youter(j)
            rdotn = rx*px0(j) + ry*py0(j)
            rdotden = rx*denx(j) + ry*deny(j)
            rho2 = rx**2.d0 + ry**2.d0

            ux(k) = ux(k) + 
     $          rdotn*rdotden/rho2/rho2*rx*speed0(j)
            uy(k) = uy(k) + 
     $          rdotn*rdotden/rho2/rho2*ry*speed0(j)
          endif
        enddo
        ux(k) = ux(k)*twopi/dble(nouter)/pi
        uy(k) = uy(k)*twopi/dble(nouter)/pi
      enddo
c     Apply trapezoid rule, but skip the diagonal term

      do k=1,nouter
        tdotden = py0(k)*den(k) - px0(k)*den(k+nouter)
        vel(k) = vel(k) + ux(k) - 
     $        cur0(k)*tdotden*py0(k)*
     $        speed0(k)/dble(nouter)
        vel(k+nouter) = vel(k + nouter) + uy(k) + 
     $        cur0(k)*tdotden*px0(k)*
     $        speed0(k)/dble(nouter)
      enddo
c     Add in correction at diagonal term that involves the curvature
c     END OF TARGET POINTS == OUTER WALL

c     START OF TARGET POINTS == OBSTACLES
      do itar = 1,nbodies
c       loop over target points
        do k=1,ninner
          ux(k) = 0.d0
          uy(k) = 0.d0

c         loop over source points
          do j=1,nouter
            rx = x((itar-1)*ninner+k) - xouter(j)
            ry = y((itar-1)*ninner+k) - youter(j)
            rdotn = rx*px0(j) + ry*py0(j)
            rdotden = rx*denx(j) + ry*deny(j)
            rho2 = rx**2.d0 + ry**2.d0

            ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*speed0(j)
            uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*speed0(j)
          enddo
          ux(k) = ux(k)*twopi/dble(nouter)/pi
          uy(k) = uy(k)*twopi/dble(nouter)/pi
        enddo
        
        do k = 1,ninner
          vel(2*nouter+(itar-1)*2*ninner+k) = 
     $      vel(2*nouter+(itar-1)*2*ninner+k) + ux(k)
          vel(2*nouter+(itar-1)*2*ninner+k+ninner) = 
     $      vel(2*nouter+(itar-1)*2*ninner+k+ninner) + uy(k)
        enddo
c       THIS SEEMS TO HAVE THE LARGEST EFFECT ON THE GMRES ITERATIONS
      enddo
c     END OF TARGET POINTS == OBSTACLE

      do k = 1,nouter
        vel(k) = vel(k) - 5.d-1*denx(k)
        vel(k + nouter) = vel(k + nouter) - 5.d-1*deny(k)
      enddo
c     add in the jump condition
c     END OF SOURCE POINTS == OUTER WALL

c************************************************************

c     START OF SOURCE POINTS == OBSTACLES
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
        enddo
c       density function due to obstacle isou

c       START OF TARGET POINTS == OBSTACLE isou
c       loop over target points
        do k=1,ninner
          ux(k) = 0.d0
          uy(k) = 0.d0

c         loop over source points
          do j=1,ninner
            if (j .ne. k) then
              rx = x((isou-1)*ninner+k) - x((isou-1)*ninner+j)
              ry = y((isou-1)*ninner+k) - y((isou-1)*ninner+j)
              rdotn = rx*px((isou-1)*ninner+j) + 
     $                ry*py((isou-1)*ninner+j)
              rdotden = rx*denx(j) + ry*deny(j)
              rho2 = rx**2.d0 + ry**2.d0

              ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*
     $            speed((isou-1)*ninner+j)
              uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*
     $            speed((isou-1)*ninner+j)
            endif
          enddo
          ux(k) = ux(k)*twopi/dble(ninner)/pi
          uy(k) = uy(k)*twopi/dble(ninner)/pi
        enddo
c       Apply trapezoid rule, but skip the diagonal term

        do k=1,ninner
          tdotden = py((isou-1)*ninner+k)*denx(k) - 
     $              px((isou-1)*ninner+k)*deny(k)
          vel(2*nouter + (isou-1)*2*ninner + k) = 
     $        vel(2*nouter + (isou-1)*2*ninner + k) +
     $        ux(k) - cur((isou-1)*ninner + k)*tdotden*
     $        py((isou-1)*ninner + k)*speed((isou-1)*ninner + k)/
     $        dble(ninner)
          vel(2*nouter + (isou-1)*2*ninner + k + ninner) = 
     $        vel(2*nouter + (isou-1)*2*ninner + k + ninner) +
     $        uy(k) + cur((isou-1)*ninner + k)*tdotden*
     $        px((isou-1)*ninner + k)*speed((isou-1)*ninner + k)/
     $        dble(ninner)
        enddo
c       Add in correction at diagonal term that involves the curvature
c       END OF TARGET POINTS == OBSTACLE isou

c       START OF TARGET POINTS ~= OBSTACLE isou
        do itar = 1,nbodies
          if (itar .eq. isou) then
            cycle
          endif
c         skip the diagonal term since this was taking care above with
c         the trapezoid rule with the correcting liming value at the
c         diagonal

c         loop over target points
          do k = 1,ninner
            ux(k) = 0.d0
            uy(k) = 0.d0

c           loop over source points
            do j = 1,ninner
              rx = x((itar-1)*ninner+k) - x((isou-1)*ninner+j)
              ry = y((itar-1)*ninner+k) - y((isou-1)*ninner+j)
              rdotn = rx*px((isou-1)*ninner+j) + 
     $                ry*py((isou-1)*ninner+j) 
              rdotden = rx*denx(j) + ry*deny(j)
              rho2 = rx**2.d0 + ry**2.d0

              ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*
     $            speed((isou-1)*ninner+j)
              uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*
     $            speed((isou-1)*ninner+j)
            enddo
            ux(k) = ux(k)*twopi/dble(ninner)/pi
            uy(k) = uy(k)*twopi/dble(ninner)/pi
          enddo
c         Apply trapezoid rule

          do k = 1,ninner
            vel(2*nouter + (itar-1)*2*ninner + k) = 
     $          vel(2*nouter + (itar-1)*2*ninner + k) + ux(k)
            vel(2*nouter + (itar-1)*2*ninner + k + ninner) = 
     $          vel(2*nouter + (itar-1)*2*ninner + k + ninner) + uy(k)
          enddo
        enddo
c       END OF TARGET POINTS ~= OBSTACLE isou


c       START OF TARGET POINTS == OUTER WALL
c       loop over target points
        do k=1,nouter
          ux(k) = 0.d0
          uy(k) = 0.d0

c         loop over source points
          do j=1,ninner
            rx = xouter(k) - x((isou-1)*ninner + j)
            ry = youter(k) - y((isou-1)*ninner + j)
            rdotn = rx*px((isou-1)*ninner + j) + 
     $              ry*py((isou-1)*ninner + j)
            rdotden = rx*denx(j) + ry*deny(j)
            rho2 = rx**2.d0 + ry**2.d0

            ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*
     $          speed((isou-1)*ninner+j)
            uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*
     $          speed((isou-1)*ninner+j)
          enddo
          ux(k) = ux(k)*twopi/dble(ninner)/pi
          uy(k) = uy(k)*twopi/dble(ninner)/pi
        enddo

        do k=1,nouter
          vel(k) = vel(k) + ux(k)
          vel(k + nouter) = vel(k + nouter) + uy(k)
        enddo
c       END OF TARGET POINTS == OUTER WALL

        do k = 1,ninner
          vel(2*nouter + (isou-1)*2*ninner + k) = 
     $    vel(2*nouter + (isou-1)*2*ninner + k) - 5.d-1*denx(k)
          vel(2*nouter + (isou-1)*2*ninner + k + ninner) = 
     $    vel(2*nouter + (isou-1)*2*ninner + k + ninner) - 5.d-1*deny(k)
        enddo
c       add in the jump condition
      enddo
c     END OF SOURCE POINTS == INNER WALL isou

c************************************************************

c     START OF SOURCE POINTS == ROTLETS AND STOKESLETS
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)

c       START OF TARGET POINTS == OUTER WALL
c       loop over target points
        do k = 1,nouter
          ux(k) = 0.d0
          uy(k) = 0.d0
          rx = xouter(k) - centerx(ibod)
          ry = youter(k) - centery(ibod)
          rho2 = rx**2.d0 + ry**2.d0
          rdots = rx*sto1 + ry*sto2
          ux(k) = 5.d-1/twopi*
     $        (-5.d-1*log(rho2)*sto1 + rdots/rho2*rx)
          uy(k) = 5.d-1/twopi*
     $        (-5.d-1*log(rho2)*sto2 + rdots/rho2*ry)
          ux(k) = ux(k) + rot*ry/rho2
          uy(k) = uy(k) - rot*rx/rho2
        enddo
c       stokeslets and rotlets contribute to the velocity

        do k=1,nouter
          vel(k) = vel(k) + ux(k)
          vel(k+nouter) = vel(k+nouter) + uy(k)
        enddo
c       END OF TARGET POINTS == OUTER WALL

c       START OF TARGET POINTS == OBSTACLE
        do itar = 1,nbodies
c         loop over target points
          do k = 1,ninner
            rx = x((itar-1)*ninner + k) - centerx(ibod)
            ry = y((itar-1)*ninner + k) - centery(ibod)
            rho2 = rx**2.d0 + ry**2.d0
            rdots = rx*sto1 + ry*sto2
            ux(k) = 5.d-1/twopi*
     $          (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
            uy(k) = 5.d-1/twopi*
     $          (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
            ux(k) = ux(k) + rot*ry/rho2
            uy(k) = uy(k) - rot*rx/rho2
          enddo
c         stokeslets and rotlets contribute to the velocity

          do k=1,ninner
            vel(2*nouter+(itar-1)*2*ninner+k) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k) + ux(k)
            vel(2*nouter+(itar-1)*2*ninner+k+ninner) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k+ninner) + uy(k)
          enddo
        enddo
c       END OF TARGET POINTS == OBSTACLE
      enddo
c     END OF SOURCE POINTS == ROTLETS AND STOKESLETS

c************************************************************

c     START OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c     STOKESLETS
      do ibod = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (ibod-1)*2*ninner + k)
          deny(k) = den(2*nouter + (ibod-1)*2*ninner + k + ninner)
        enddo
c       density function due to obstacle ibod

        do k = 1,ninner
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) - 
     $      denx(k)*speed((ibod-1)*ninner+k)
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) - 
     $      deny(k)*speed((ibod-1)*ninner+k)
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) - 
     $      (denx(k)*y((ibod-1)*ninner+k)-
     $      deny(k)*x((ibod-1)*ninner+k))*
     $      speed((ibod-1)*ninner+k)
        enddo
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1)/twopi * 
     $    (twopi/dble(ninner))
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2)/twopi * 
     $    (twopi/dble(ninner))
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3)/twopi * 
     $    (twopi/dble(ninner))

        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) + sto1
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) + sto2
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) + rot
c       END OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c       STOKESLETS
      enddo

c************************************************************

c     START OF RANK 1 CORRECTION TO HANDLE NULL SPACE FROM THE OUTERMOST
c     BOUNDARY
      do k = 1,nouter
        denx(k) = den(k)
        deny(k) = den(k + nouter)
      enddo

      sigdotn = 0.d0
      do k = 1,nouter
        sigdotn = sigdotn + (px0(k)*denx(k) + py0(k)*deny(k))*speed0(k)
      enddo
      sigdotn = sigdotn*twopi/dble(nouter)
c     integral of dot product of the normal vector and density function
c     defined on the outer boundary

      do k = 1,nouter
        vel(k) = vel(k) + sigdotn * px0(k)
        vel(k+nouter) = vel(k+nouter) + sigdotn * py0(k)
      enddo
c     outer product of normal at source with normal at target multiplied
c     by density function.  This removes the rank one null space


      do k = 1,2*nouter + 2*ninner*nbodies
        print*,vel(k)
      enddo

      return
      end
c***********************************************************************
      subroutine matvec_DLP_fmm(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 50)

      dimension den(ntotal)
      dimension vel(ntotal)

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter

      dimension x(nmax),y(nmax)
      dimension px(nmax),py(nmax)
      dimension cur(nmax),speed(nmax)
      dimension centerx(maxbodies),centery(maxbodies)

      dimension xouter(nmax),youter(nmax)
      dimension px0(nmax),py0(nmax)
      dimension cur0(nmax),speed0(nmax)

      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))
      dimension ux(max(ninner,nouter)),uy(max(ninner,nouter))

      real *8 xboth(nouter + ninner*nbodies)
      real *8 yboth(nouter + ninner*nbodies)
      complex *16 eye
      complex *16 mu(nouter + ninner*nbodies)
      complex *16 z(nouter + ninner*nbodies)
      complex *16 dip1(nouter + ninner*nbodies)
      complex *16 dip2(nouter + ninner*nbodies)
      complex *16 vel_cmplx(nouter + ninner*nbodies)

      czero = (0.0d0,0.0d0)
      eye = (0.d0,1.d0)
      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do k = 1,nouter + ninner*nbodies
        vel_cmplx(k) = czero
      enddo
c     initialize complex-valued velocity to 0

      do k = 1,2*nouter + 2*ninner*nbodies + 3*nbodies
        vel(k) = 0.d0
      enddo
c     initialize velocity to 0

      do k = 1,nouter
        mu(k) = (den(k+nouter) - eye*den(k))*speed0(k)*
     $      twopi/dble(nouter)
        z(k) = px0(k) + eye*py0(k)
        xboth(k) = xouter(k)
        yboth(k) = youter(k)
      enddo

      do isou = 1,nbodies
        do k=1,ninner
          mu(nouter + (isou-1)*ninner + k) = 
     $        (den(2*nouter+(isou-1)*2*ninner+k+ninner) - eye*
     $         den(2*nouter+(isou-1)*2*ninner+k))*
     $         speed((isou-1)*ninner+k)*twopi/dble(ninner)
          z(nouter + (isou-1)*ninner + k) = 
     $        px((isou-1)*ninner+k) + eye*py((isou-1)*ninner+k)
          xboth(nouter + (isou-1)*ninner + k) = 
     $        x((isou-1)*ninner + k)
          yboth(nouter + (isou-1)*ninner + k) = 
     $        y((isou-1)*ninner + k)
        enddo
      enddo

      dip1 = 2.5d-1/pi*mu*z
      dip2 = 2.5d-1/pi*(mu*conjg(z) - conjg(mu)*z)

      call stokesDLPnew(nouter+nbodies*ninner,xboth,yboth,
     $        dip1,dip2,vel_cmplx)
c     Apply FMM to do the all-to-all particle interactions in linear
c     time

      do k = 1,nouter
        vel(k) = -imag(vel_cmplx(k))
        vel(k+nouter) = real(vel_cmplx(k))
      enddo
c     copy velocity on outer boundary in correct order

      do isou = 1,nbodies
        do k = 1,ninner
          vel(2*nouter+(isou-1)*2*ninner+k) = 
     $      -imag(vel_cmplx(nouter+(isou-1)*ninner+k))
          vel(2*nouter+(isou-1)*2*ninner+k+ninner) = 
     $      real(vel_cmplx(nouter+(isou-1)*ninner+k))
        enddo
      enddo
c     copy velocity on inner boundaries in correct order


c     START OF EXTRA TERMS DUE TO CHARGES ON OUTER WALLS
      do k = 1,nouter
        denx(k) = den(k)
        deny(k) = den(k + nouter)
      enddo

      do k=1,nouter
        tdotden = py0(k)*denx(k) - px0(k)*deny(k)
        vel(k) = vel(k) - 
     $        cur0(k)*tdotden*py0(k)*
     $        speed0(k)/dble(nouter)
        vel(k+nouter) = vel(k + nouter) + 
     $        cur0(k)*tdotden*px0(k)*
     $        speed0(k)/dble(nouter)
      enddo
c     Add in correction at diagonal term that involves the curvature

      do k = 1,nouter
        vel(k) = vel(k) - 5.d-1*denx(k)
        vel(k + nouter) = vel(k + nouter) - 5.d-1*deny(k)
      enddo
c     add in the jump condition
c     END OF EXTRA TERMS DUE TO CHARGES ON OUTER WALLS

c************************************************************

c     START OF EXTRA TERMS DUE TO CHARGES ON OBSTACLE isou
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
        enddo

        do k=1,ninner
          tdotden = py((isou-1)*ninner+k)*denx(k) - 
     $              px((isou-1)*ninner+k)*deny(k)
          vel(2*nouter + (isou-1)*2*ninner + k) = 
     $        vel(2*nouter + (isou-1)*2*ninner + k) -
     $        cur((isou-1)*ninner + k)*tdotden*
     $        py((isou-1)*ninner + k)*speed((isou-1)*ninner + k)/
     $        dble(ninner)
          vel(2*nouter + (isou-1)*2*ninner + k + ninner) = 
     $        vel(2*nouter + (isou-1)*2*ninner + k + ninner) +
     $        cur((isou-1)*ninner + k)*tdotden*
     $        px((isou-1)*ninner + k)*speed((isou-1)*ninner + k)/
     $        dble(ninner)
        enddo
c       Add in correction at diagonal term that involves the curvature

        do k = 1,ninner
          vel(2*nouter + (isou-1)*2*ninner + k) = 
     $    vel(2*nouter + (isou-1)*2*ninner + k) - 5.d-1*denx(k)
          vel(2*nouter + (isou-1)*2*ninner + k + ninner) = 
     $    vel(2*nouter + (isou-1)*2*ninner + k + ninner) - 5.d-1*deny(k)
        enddo
c       add in the jump condition
      enddo
c     END OF EXTRA TERMS DUE TO CHARGES ON OBSTACLE isou



c************************************************************
c     START OF SOURCE POINTS == ROTLETS AND STOKESLETS
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)

c       START OF TARGET POINTS == OUTER WALL
c       loop over target points
        do k = 1,nouter
          ux(k) = 0.d0
          uy(k) = 0.d0
          rx = xouter(k) - centerx(ibod)
          ry = youter(k) - centery(ibod)
          rho2 = rx**2.d0 + ry**2.d0
          rdots = rx*sto1 + ry*sto2
          ux(k) = 5.d-1/twopi*
     $        (-5.d-1*log(rho2)*sto1 + rdots/rho2*rx)
          uy(k) = 5.d-1/twopi*
     $        (-5.d-1*log(rho2)*sto2 + rdots/rho2*ry)
          ux(k) = ux(k) + rot*ry/rho2
          uy(k) = uy(k) - rot*rx/rho2
        enddo
c       stokeslets and rotlets contribute to the velocity

        do k=1,nouter
          vel(k) = vel(k) + ux(k)
          vel(k+nouter) = vel(k+nouter) + uy(k)
        enddo
c       END OF TARGET POINTS == OUTER WALL

c       START OF TARGET POINTS == OBSTACLE
        do itar = 1,nbodies
c         loop over target points
          do k = 1,ninner
            ux(k) = 0.d0
            uy(k) = 0.d0
            rx = x((itar-1)*ninner + k) - centerx(ibod)
            ry = y((itar-1)*ninner + k) - centery(ibod)
            rho2 = rx**2.d0 + ry**2.d0
            rdots = rx*sto1 + ry*sto2
            ux(k) = 5.d-1/twopi*
     $          (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
            uy(k) = 5.d-1/twopi*
     $          (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
            ux(k) = ux(k) + rot*ry/rho2
            uy(k) = uy(k) - rot*rx/rho2
          enddo
c         stokeslets and rotlets contribute to the velocity

          do k=1,ninner
            vel(2*nouter+(itar-1)*2*ninner+k) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k) + ux(k)
            vel(2*nouter+(itar-1)*2*ninner+k+ninner) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k+ninner) + uy(k)
          enddo
        enddo
c       END OF TARGET POINTS == OBSTACLE
      enddo
c     END OF SOURCE POINTS == ROTLETS AND STOKESLETS

c************************************************************

c     START OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c     STOKESLETS
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)
        do k = 1,ninner
          denx(k) = den(2*nouter + (ibod-1)*2*ninner + k)
          deny(k) = den(2*nouter + (ibod-1)*2*ninner + k + ninner)
        enddo
c       density function due to obstacle ibod

        do k = 1,ninner
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) - 
     $      denx(k)*speed((ibod-1)*ninner+k)
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) - 
     $      deny(k)*speed((ibod-1)*ninner+k)
          vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) - 
     $      (denx(k)*y((ibod-1)*ninner+k)-
     $      deny(k)*x((ibod-1)*ninner+k))*
     $      speed((ibod-1)*ninner+k)
        enddo
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1)/twopi * 
     $    (twopi/dble(ninner))
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2)/twopi * 
     $    (twopi/dble(ninner))
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $    vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3)/twopi * 
     $    (twopi/dble(ninner))

        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+1) + sto1
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+2) + sto2
        vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) = 
     $      vel(2*nouter+2*nbodies*ninner+(ibod-1)*3+3) + rot
c       END OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c       STOKESLETS
      enddo

c************************************************************

c     START OF RANK 1 CORRECTION TO HANDLE NULL SPACE FROM THE OUTERMOST
c     BOUNDARY
      do k = 1,nouter
        denx(k) = den(k)
        deny(k) = den(k + nouter)
      enddo

      sigdotn = 0.d0
      do k = 1,nouter
        sigdotn = sigdotn + (px0(k)*denx(k) + py0(k)*deny(k))*speed0(k)
      enddo
      sigdotn = sigdotn*twopi/dble(nouter)
c     integral of dot product of the normal vector and density function
c     defined on the outer boundary

      do k = 1,nouter
        vel(k) = vel(k) + sigdotn * px0(k)
        vel(k+nouter) = vel(k+nouter) + sigdotn * py0(k)
      enddo
c     outer product of normal at source with normal at target multiplied
c     by density function.  This removes the rank one null space


      return
      end

c***********************************************************************
      subroutine msolve_DLP(nn,r,z,nelt,ia,ja,a,isym,rwork,iwork)
c     Can put preconditioner in this routine.  For now, use the identity
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)

      dimension r(nn),z(nn)
      dimension iwork(3)

      complex *16 eye
      real *8 wsave(4*nmax+15)
      complex *16 zden1(nmax),zden2(nmax)
      complex *16 g_minus1,g_plus1
      complex *16 h_minus1,h_plus1
      complex *16 alpha_minus1,alpha_plus1
      complex *16 beta_minus1,beta_plus1

      eye = (0.0d0,1.0d0)

      ninner = iwork(1)
      nbodies = iwork(2)
      nouter = iwork(3)

      do i = 1,2*nouter
        z(i) = r(i)
      enddo
c     identity preconditioner for outer boundary

      do i = 2*nouter+2*nbodies*ninner+1,nn
        z(i) = r(i)
      enddo
c     identity preconditioner for rotlet and stokeslet terms

c     START APPLYING PRECONDITIONER TO INNER BOUNDARY
      call DCFFTI(ninner,wsave)
      do k = 1,nbodies
        do j = 1,ninner
          zden1(j) = r(2*nouter + (k-1)*2*ninner + j)
          zden2(j) = r(2*nouter + (k-1)*2*ninner + j + ninner)
        enddo
        call DCFFTF(ninner,zden1,wsave)
        call DCFFTF(ninner,zden2,wsave)
c       computer Fourier coefficients of the density function on
c       the current body.  Need to divide by number of points since
c       DCFFTF and DCFFTB are not inverses, but scale up by the number
c       of points
        zden1 = zden1/dble(ninner)
        zden2 = zden2/dble(ninner)
        g_minus1 = zden1(ninner)
        h_minus1 = zden2(ninner)
        g_plus1 = zden1(2)
        h_plus1 = zden2(2)

        do j=2,ninner
          zden1(j) = -2.d0*zden1(j)
          zden2(j) = -2.d0*zden2(j)
        enddo
c       invert the Fourier modies with frequency greater than or equal
c       to 2.  The DLP of these modes is 0, so only the -1/2*identity
c       term plays a role

        alpha_minus1 = -1.25d0*g_minus1 - 7.5d-1*eye*h_minus1  
     $                 -2.5d-1*g_plus1 - 2.5d-1*eye*h_plus1
        beta_minus1 = +7.5d-1*eye*g_minus1 - 1.25d0*h_minus1  
     $                -2.5d-1*eye*g_plus1 + 2.5d-1*h_plus1
        alpha_plus1 = -2.5d-1*g_minus1 + 2.5d-1*eye*h_minus1  
     $                -1.25d0*g_plus1 + 7.5d-1*eye*h_plus1
        beta_plus1 = +2.5d-1*eye*g_minus1 + 2.5d-1*h_minus1  
     $               -7.5d-1*eye*g_plus1 - 1.25d0*h_plus1
c       Use pseudo-inverse of the lienar system that relates the one and
c       minus one nodes of the density function to the one and minus one
c       nodes of the right hand side

        zden1(2) = alpha_plus1
        zden1(ninner) = alpha_minus1
        zden2(2) = beta_plus1
        zden2(ninner) = beta_minus1
c       Assign computed Fourier nodes to appropriate locations of zden1
c       and zden2

        call DCFFTB(ninner,zden1,wsave)
        call DCFFTB(ninner,zden2,wsave)
c       move back to physical space
        do j = 1,ninner
          z(2*nouter + (k-1)*2*ninner + j) = real(zden1(j))
          z(2*nouter + (k-1)*2*ninner + j + ninner) = real(zden2(j))
        enddo
c       assign real and imaginary parts to appropriate locations of r
      enddo

      return
      end


c***********************************************************************
c     TODO: NEW INPUT VARIABLE FOR THE BARYCENTRIC FLAG
      subroutine deformation_on_boundary(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor at points on the boundary of the
c     obstacle
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xouter(nouter),youter(nouter)
      dimension px0(nouter),py0(nouter)
      dimension speed0(nouter)
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)

      dimension E11(ninner*nbodies)
      dimension E12(ninner*nbodies)
      dimension E22(ninner*nbodies)

      complex *16 zden(ninner)
      real*8 wsave(4*ninner+15)
      complex *16 eye
      real *8 denx(max(ninner,nouter)),deny(max(ninner,nouter))

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      eye = (0.d0,1.d0)

      do j=1,ninner*nbodies
        E11(j) = 0.d0
        E12(j) = 0.d0
        E22(j) = 0.d0
      enddo
c     initialize components to be zero


c     START OF SOURCE POINTS ==  SOLID WALL, TARGET POINTS == OBSTACLES
      do j=1,nouter
        denx(j) = den(j)
        deny(j) = den(j+nouter)
      enddo
c     Density function defined on the outer geometry

      do itar=1,nbodies
c       loop over targets
        do k=1,ninner
c         loop over sources
          do j=1,nouter
            rx = x((itar-1)*ninner+k) - xouter(j)
            ry = y((itar-1)*ninner+k) - youter(j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px0(j) + ry*py0(j)
            routn = rx*py0(j) + ry*px0(j)
            rdotden = rx*denx(j) + ry*deny(j)
            routden = rx*deny(j) + ry*denx(j)

            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + 5.d-1*
     $          (2.d0*rdotn*rdotden/rho2/rho2 + 
     $          2.d0*rdotden/rho2/rho2*rx*px0(j) +
     $          2.d0*rdotn/rho2/rho2*rx*denx(j) -
     $          8.d0*rdotn*rdotden/rho2**3.d0*rx*rx)*
     $          speed0(j)*twopi/dble(nouter)/pi

            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 5.d-1*
     $          (rdotden*routn/rho2/rho2 +
     $          rdotn*routden/rho2/rho2 -
     $          8.d0*rdotn*rdotden*rx*ry/rho2**3.d0)*
     $          speed0(j)*twopi/dble(nouter)/pi

            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 5.d-1*
     $          (2.d0*rdotn*rdotden/rho2/rho2 + 
     $          2.d0*rdotden/rho2/rho2*ry*py0(j) +
     $          2.d0*rdotn/rho2/rho2*ry*deny(j) -
     $          8.d0*rdotn*rdotden/rho2**3.d0*ry*ry)*
     $          speed0(j)*twopi/dble(nouter)/pi
          enddo
        enddo
      enddo
c     END OF SOURCE POINTS ==  SOLID WALL, TARGET POINTS == OBSTACLES

c     START OF SOURCE POINTS == OBSTACLES, TARGET POINTS == OBSTACLE
c     isou
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
        enddo
c       density function due to obstacle isou

c       START OF TARGET POINTS == OBSTACLE isou
c       loop over target points
        do k = 1,ninner,2
c         loop over source points
          do j = 2,ninner,2
            rx = x((isou-1)*ninner+k) - x((isou-1)*ninner+j)
            ry = y((isou-1)*ninner+k) - y((isou-1)*ninner+j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner+j) + 
     $              ry*py((isou-1)*ninner+j)
            routn = rx*py((isou-1)*ninner+j) + 
     $              ry*px((isou-1)*ninner+j)

            sx = denx(j) - denx(k)
            sy = deny(j) - deny(k)
c           subtract of the density function evaluated at the target point
c           to reduce the singularity to 1/r.  Then, odd-even integration
c           converges to the PV integral
            rdots = rx*sx + ry*sy
            routs = rx*sy + ry*sx

            E11((isou-1)*ninner+k) = E11((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (2.d0*rdotn*rdots/rho2/rho2 + 
     $          2.d0*rdots/rho2/rho2*rx*px((isou-1)*ninner+j) +
     $          2.d0*rdotn/rho2/rho2*rx*sx -
     $          8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

            E12((isou-1)*ninner+k) = E12((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (rdots*routn/rho2/rho2 +
     $          rdotn*routs/rho2/rho2 -
     $          8.d0*rdotn*rdots*rx*ry/rho2**3.d0)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

            E22((isou-1)*ninner+k) = E22((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (2.d0*rdotn*rdots/rho2/rho2 + 
     $          2.d0*rdots/rho2/rho2*ry*py((isou-1)*ninner+j) +
     $          2.d0*rdotn/rho2/rho2*ry*sy -
     $          8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
          enddo
        enddo
c       need to multiply by 2 since the grid spacing is twice as large
c       compute deformation tensor on odd indexed terms

c       loop over target points
        do k = 2,ninner,2
c         loop over source points
          do j = 1,ninner,2
            rx = x((isou-1)*ninner+k) - x((isou-1)*ninner+j)
            ry = y((isou-1)*ninner+k) - y((isou-1)*ninner+j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner+j) + 
     $              ry*py((isou-1)*ninner+j)
            routn = rx*py((isou-1)*ninner+j) + 
     $              ry*px((isou-1)*ninner+j)

            sx = denx(j) - denx(k)
            sy = deny(j) - deny(k)
c           subtract of the density function evaluated at the target point
c           to reduce the singularity to 1/r.  Then, odd-even integration
c           converges to the PV integral
            rdots = rx*sx + ry*sy
            routs = rx*sy + ry*sx

            E11((isou-1)*ninner+k) = E11((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (2.d0*rdotn*rdots/rho2/rho2 + 
     $          2.d0*rdots/rho2/rho2*rx*px((isou-1)*ninner+j) +
     $          2.d0*rdotn/rho2/rho2*rx*sx -
     $          8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

            E12((isou-1)*ninner+k) = E12((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (rdots*routn/rho2/rho2 +
     $          rdotn*routs/rho2/rho2 -
     $          8.d0*rdotn*rdots*rx*ry/rho2**3.d0)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

            E22((isou-1)*ninner+k) = E22((isou-1)*ninner+k) + 
     $          2.d0*5.d-1*
     $          (2.d0*rdotn*rdots/rho2/rho2 + 
     $          2.d0*rdots/rho2/rho2*ry*py((isou-1)*ninner+j) +
     $          2.d0*rdotn/rho2/rho2*ry*sy -
     $          8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
     $          speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
          enddo
        enddo
c       need to multiply by 2 since the grid spacing is twice as large
c       compute deformation tensor on odd indexed terms
c       END OF SOURCE POINTS == OBSTACLES, TARGET POINTS == OBSTACLE
c       isou


c       START OF SOURCE POINTS == OBSTACLES, TARGET POINTS ~= OBSTACLE
c       isou
        do itar = 1,nbodies
          if (itar .eq. isou) then
            cycle
          endif
c         skip the diagonal term since this was taking care above with
c         the trapezoid rule with odd even integration

c         loop over target points
          do k = 1,ninner
c           loop over source points
            do j = 1,ninner
              rx = x((itar-1)*ninner+k) - x((isou-1)*ninner+j)
              ry = y((itar-1)*ninner+k) - y((isou-1)*ninner+j)
              rho2 = rx**2.d0 + ry**2.d0
              rdotn = rx*px((isou-1)*ninner+j) +
     $                ry*py((isou-1)*ninner+j)
              routn = rx*py((isou-1)*ninner+j) +
     $                ry*px((isou-1)*ninner+j)
              sx = denx(j)
              sy = deny(j)
              rdots = rx*sx + ry*sy
              routs = rx*sy + ry*sx

              E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + 
     $            5.d-1*
     $            (2.d0*rdotn*rdots/rho2/rho2 + 
     $            2.d0*rdots/rho2/rho2*rx*px((isou-1)*ninner+j) +
     $            2.d0*rdotn/rho2/rho2*rx*sx -
     $            8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

              E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 
     $            5.d-1*
     $            (rdots*routn/rho2/rho2 +
     $            rdotn*routs/rho2/rho2 -
     $            8.d0*rdotn*rdots*rx*ry/rho2**3.d0)*
     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi

              E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 
     $            5.d-1*
     $            (2.d0*rdotn*rdots/rho2/rho2 + 
     $            2.d0*rdots/rho2/rho2*ry*py((isou-1)*ninner+j) +
     $            2.d0*rdotn/rho2/rho2*ry*sy -
     $            8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
            enddo
          enddo
        enddo
c       END OF SOURCE POINTS == OBSTACLES, TARGET POINTS ~= OBSTACLE
c       isou
      enddo

c     START OF JUMP ALONG DIAGONAL TERM
      call DCFFTI(ninner,wsave)
      do ibod = 1,nbodies
        do k=1,ninner 
          zden(k) = den(2*nouter+(ibod-1)*2*ninner+k) + 
     $          eye*den(2*nouter+(ibod-1)*2*ninner+k+ninner)
        enddo
        call fourierDiff(ninner,zden,wsave)
c       real part of zden is parameter derivative of the first component
c       of the density function and the complex part is the derivative of
c       the second component of the density function

        do k=1,ninner
          tx = -py((ibod-1)*ninner+k)
          ty =  px((ibod-1)*ninner+k)
          dsdtx = dreal(zden(k))/speed((ibod-1)*ninner+k)
          dsdty = dimag(zden(k))/speed((ibod-1)*ninner+k)
          dsdt_dot_tau = dsdtx*tx + dsdty*ty
          E11((ibod-1)*ninner+k) = E11((ibod-1)*ninner+k) + 5.d-1*
     $        dsdt_dot_tau*(tx**2.d0 - ty**2.d0)
          E12((ibod-1)*ninner+k) = E12((ibod-1)*ninner+k) + 5.d-1*
     $        dsdt_dot_tau*2.d0*tx*ty
          E22((ibod-1)*ninner+k) = E22((ibod-1)*ninner+k) + 5.d-1*
     $        dsdt_dot_tau*(-tx**2.d0 + ty**2.d0)
        enddo
c       Add in jump conditions
      enddo
c     END OF JUMP ALONG DIAGONAL TERM


c     START OF SOURCE POINTS == ROTLETS AND STOKESLETS, TARGET POINTS ==
c     OBSTACLE
      do ibod = 1,nbodies
        sto1 = den(2*nouter + 2*ninner*nbodies + (ibod-1)*3+1)
        sto2 = den(2*nouter + 2*ninner*nbodies + (ibod-1)*3+2)
        rot  = den(2*nouter + 2*ninner*nbodies + (ibod-1)*3+3)
c       loop over targets
        do itar = 1,nbodies
c        do itar = ibod,ibod
          do k = 1,ninner
            rx = x((itar-1)*ninner+k) - centerx(ibod)
            ry = y((itar-1)*ninner+k) - centery(ibod)
            rho2 = rx**2.d0 + ry**2.d0
            rdots = rx*sto1 + ry*sto2
            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + 
     $          5.d-1/twopi*(rdots/rho2 - 2*rx*rx*rdots/rho2/rho2)
            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) - 
     $          1.d0/twopi*(rdots/rho2/rho2*rx*ry)
            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 
     $          5.d-1/twopi*(rdots/rho2 - 2*ry*ry*rdots/rho2/rho2)
c         stokeslet contribution

            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) - 
     $          2.d0*rot*rx*ry/rho2/rho2
            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 
     $          rot*(rx**2.d0 - ry**2.d0)/rho2/rho2
            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) - 
     $          2.d0*rot*rx*ry/rho2/rho2
c         rotlet contribution
          enddo
        enddo
c       Add in contribution from rotlets and stokeslets
      enddo
c     END OF SOURCE POINTS == ROTLETS AND STOKESLETS, TARGETS ==
c     OBSTACLES

      end


c***********************************************************************
      subroutine computeShearStress(ninner,nbodies,nouter,x,y,den,
     $    shear_stress)
c     Input x and y coordinates and the density function on the
c     boundary and return the shear stress on the inner walls.  
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
c     x and y coordinates of the normal of the obstacle
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
c     Jacobian and curvature of the geometry

      dimension xouter(nouter*nbodies),youter(nouter*nbodies)
c     x and y coordinates of confining wall
      dimension px0(nouter*nbodies), py0(nouter*nbodies)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nouter*nbodies), speed0(nouter*nbodies)
c     Jacobian and curvature of the confining wall

      dimension den(2*ninner*nbodies + 3*nbodies + 2*nouter)
c     x and y coordinates of the density function

      dimension E11(ninner*nbodies)
      dimension E12(ninner*nbodies)
      dimension E22(ninner*nbodies)
c     components of the deformation gradient on the obstacle
      dimension tractionx(ninner),tractiony(ninner)
c     traction which is computed one body at a time
      dimension shear_stress(ninner*nbodies)
c     shear stress on the obstacle

c     can't declare input variables into a common field, so need new
c     variable name for x,y,ninner,nbodies
      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c
      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     load geometry of initial shape

      call deformation_on_boundary(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor on the boundary of the interface


      do ibod = 1,nbodies
        do k = 1,ninner
          tractionx(k) = -2.d0*E11((ibod-1)*ninner+k)*
     $       px((ibod-1)*ninner+k) - 
     $                    2.d0*E12((ibod-1)*ninner+k)*
     $       py((ibod-1)*ninner+k)
          tractiony(k) = -2.d0*E12((ibod-1)*ninner+k)*
     $       px((ibod-1)*ninner+k) - 
     $                    2.d0*E22((ibod-1)*ninner+k)*
     $       py((ibod-1)*ninner+k)
        enddo

        do k = 1,ninner
          shear_stress((ibod-1)*ninner+k) = 
     $      tractionx(k)*py((ibod-1)*ninner+k) - 
     $      tractiony(k)*px((ibod-1)*ninner+k) 
        enddo
      enddo


      end

c***********************************************************************
c      subroutine computeQoiTargetsOld(ninner,nbodies,nouter,
c     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
cc     Compute the velocity, pressure, and vorticity at a set of target 
cc     points xtar and ytar.  Target points must be sufficiently far 
cc     away from the each boundary since no near-singular integration 
cc     is used.  
cc     It is just the vanilla trapezoid rule, but will assign a value of
cc     zero for points that are inside or close to the boundary
c      implicit real*8 (a-h,o-z)
c
cc     input variables
c      dimension x(ninner*nbodies),y(ninner*nbodies)
c      dimension px(ninner*nbodies),py(ninner*nbodies)
c      dimension cur(ninner*nbodies),speed(ninner*nbodies)
c      dimension centerx(nbodies),centery(nbodies)
c      dimension xouter(nouter),youter(nouter)
cc     x and y coordinates of confining wall
c      dimension px0(nouter), py0(nouter)
cc     x and y coordinates of the normal of the confining wall
c      dimension cur0(nouter), speed0(nouter)
cc     Jacobian and curvature of the confining wall
c      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
c      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))
c
cc     output variables
c      dimension xtar(ntargets),ytar(ntargets)
c      dimension utar(ntargets),vtar(ntargets)
c      dimension press_tar(ntargets)
c      dimension vort_tar(ntargets)
c
cc     local variables
c      dimension iside(ntargets)
c      dimension nnear(nbodies)
c      dimension indnear(ntargets,nbodies)
c      dimension iup(ntargets,nbodies)
c
cc     nnear is the number of near points for each body
cc     indnear is the index of the points that are close to a 
cc     particular body
c      call classifyPoints(ninner,nbodies,x,y,
c     $      ntargets,xtar,ytar,iside,nnear,indnear,iup)
c
c      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
c     $    centerx,centery)
cc     build inner geometry
c
c      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
cc     build outer geometry
c
c      pi = 4.d0*datan(1.d0)
c      twopi = 2.d0*pi
c
c      do j = 1,ntargets
c        utar(j) = 0.d0
c        vtar(j) = 0.d0
c        press_tar(j) = 0.d0
c        vort_tar(j) = 0.d0
c      enddo
c
cc     Contribution from the density function on the solid wall
c      do k=1,nouter
c        denx(k) = den(k)
c        deny(k) = den(k+nouter)
c      enddo
cc     Density function defined on the outer geometry
c
c      do j = 1,ntargets
c        do k = 1,nouter
c          rx = xtar(j) - xouter(k)
c          ry = ytar(j) - youter(k)
c          rho2 = rx**2.d0 + ry**2.d0
c          rdotn = rx*px0(k) + ry*py0(k)
c          rdotden = rx*denx(k) + ry*deny(k)
c          dendotn = px0(k)*denx(k) + py0(k)*deny(k)
c          routn = -rx*py0(k) + ry*px0(k)
c          routden = -rx*deny(k) + ry*denx(k)
c
c          utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
c     $      speed0(k)*twopi/dble(nouter)/pi
c          vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
c     $      speed0(k)*twopi/dble(nouter)/pi
c          press_tar(j) = press_tar(j) + 
c     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
c     $        speed0(k)*twopi/dble(nouter)/pi
c          vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
c     $      rho2**2.d0*speed0(k)*twopi/dble(nouter)/pi
c        enddo
c      enddo
c
cc     Contribution from the density function on each inner obstacle
c      do isou = 1,nbodies
c        do k = 1,ninner
c          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
c          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
c        enddo
cc       Density function defined on the inner geometry
c
c        do j = 1,ntargets
c          do k = 1,ninner
c            rx = xtar(j) - x((isou-1)*ninner + k)
c            ry = ytar(j) - y((isou-1)*ninner + k)
c            rho2 = rx**2.d0 + ry**2.d0
c            rdotn = rx*px((isou-1)*ninner + k) + 
c     $              ry*py((isou-1)*ninner + k)
c            rdotden = rx*denx(k) + ry*deny(k)
c            dendotn = px((isou-1)*ninner + k)*denx(k) + 
c     $                py((isou-1)*ninner + k)*deny(k)
c            routn = -rx*py((isou-1)*ninner+k) + ry*px((isou-1)*ninner+k)
c            routden = -rx*deny(k) + ry*denx(k)
c
c            utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
c     $        speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
c            vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
c     $        speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
c            press_tar(j) = press_tar(j) + 
c     $          (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
c     $          speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
c            vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
c     $        rho2**2.d0*speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi 
c          enddo
c        enddo
c      enddo
c
cc     Contribution from Rotlets and Stokeslets
c      do ibod = 1,nbodies
c        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
c        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
c        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)
c        do j = 1,ntargets
c          rx = xtar(j) - centerx(ibod)
c          ry = ytar(j) - centery(ibod)
c          rho2 = rx**2.d0 + ry**2.d0
c          rdots = rx*sto1 + ry*sto2
c          utar(j) = utar(j) + 5.d-1/twopi*
c     $        (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
c          vtar(j) = vtar(j) + 5.d-1/twopi*
c     $        (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
c          press_tar(j) = press_tar(j) + 1.d0/twopi*
c     $        rdots/rho2
c          vort_tar(j) = vort_tar(j) + (ry*sto1 - rx*sto2)/
c     $        rho2/twopi
cc         stokeslet contribution
c
c          utar(j) = utar(j) + rot*ry/rho2
c          vtar(j) = vtar(j) - rot*rx/rho2
cc         rotlet contribution
c        enddo
c      enddo
cc     Add in contribution from Rotlets and Stokeslets
c
c      
c      
c
c      do k = 1,ntargets
cc        if (iside(k) .eq. 0 .or. inear(k) .eq. 1) then
c        if (iside(k) .eq. 0) then
c          utar(k) = 0.d0
c          vtar(k) = 0.d0
c          press_tar(k) = 0.d0
c          vort_tar(k) = 0.d0
c        endif
c      enddo
c
c
c      do j=1,nbodies
c        do k = 1,nnear(j)
c          if (iside(indnear(k,j)) .ne. 0) then
c            utar(indnear(k,j)) = dble(iup(k,j))
c            vtar(indnear(k,j)) = dble(iup(k,j))
c          endif
c        enddo
c      enddo
c
c
c
c      end

c***********************************************************************
c     TODO: NEW INPUT VARIABLE FOR THE BARYCENTRIC FLAG
      subroutine computeQoiTargets(ninner,nbodies,nouter,
     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
c     Compute the velocity, pressure, and vorticity at a set of target 
c     points xtar and ytar.
      implicit real*8 (a-h,o-z)

c     input variables
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xouter(nouter),youter(nouter)
c     x and y coordinates of confining wall
      dimension px0(nouter), py0(nouter)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nouter), speed0(nouter)
c     Jacobian and curvature of the confining wall
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))

c     output variables
      dimension xtar(ntargets),ytar(ntargets)
      dimension utar(ntargets),vtar(ntargets)
      dimension press_tar(ntargets)
      dimension vort_tar(ntargets)

c     local variables
      dimension iside(ntargets)
      dimension iup(ntargets)
      dimension xtar_level(ntargets),ytar_level(ntargets)
      dimension utar_level(ntargets),vtar_level(ntargets)
      dimension press_tar_level(ntargets),vort_tar_level(ntargets)
      dimension xup(16*max(ninner,nouter))
      dimension yup(16*max(ninner,nouter))
      dimension pxup(16*max(ninner,nouter))
      dimension pyup(16*max(ninner,nouter))
      dimension curup(16*max(ninner,nouter))
      dimension speedup(16*max(ninner,nouter))
      dimension denxup(16*max(ninner,nouter))
      dimension denyup(16*max(ninner,nouter))
      complex *16 eye
      complex *16 zin(16*max(ninner,nouter))
      complex *16 zout(16*max(ninner,nouter))

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      eye = (0.d0,1.d0)

      do j = 1,ntargets
        utar(j) = 0.d0
        vtar(j) = 0.d0
        press_tar(j) = 0.d0
        vort_tar(j) = 0.d0
      enddo

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

c     build inner geometry
      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     build outer geometry

c     Contribution from the density function on the solid wall
      do k=1,nouter
        denx(k) = den(k)
        deny(k) = den(k+nouter)
      enddo
c     Density function defined on the outer geometry


      do j = 1,nouter
        xup(j) = xouter(j)
        yup(j) = youter(j)
        denxup(j) = denx(j)
        denyup(j) = deny(j)
      enddo

      call classifyPoints(nouter,xup,yup,
     $    ntargets,xtar,ytar,iside,nnear,iup)

      do k = 1,ntargets
        if (iside(k) .eq. 1 .or. iup(k) .gt. 16) then
          utar(k) = 1.0d20
          vtar(k) = 1.0d20
          press_tar(k) = 1.0d20
          vort_tar(k) = 1.0d20
        endif
      enddo

c     START OF DOING POINTS THAT REQUIRE NO UPSAMPLING
      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 1 .and. iside(j) .eq. 0) then
          xtar_level(icount) = xtar(j)
          ytar_level(icount) = ytar(j)
          icount = icount + 1
        endif
      enddo

      call evalLPouter(nouter,xup,yup,
     $  denxup,denyup,icount-1,xtar_level,ytar_level,
     $  utar_level,vtar_level,press_tar_level,vort_tar_level)

      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 1 .and. iside(j) .eq. 0) then
          utar(j) = utar(j) + utar_level(icount)
          vtar(j) = vtar(j) + vtar_level(icount)
          press_tar(j) = press_tar(j) + press_tar_level(icount)
          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
          icount = icount + 1
        endif
      enddo
c     END OF DOING POINTS THAT REQUIRE NO UPSAMPLING

c     START OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 2 .and. iside(j) .eq. 0) then
          xtar_level(icount) = xtar(j)
          ytar_level(icount) = ytar(j)
          icount = icount + 1
        endif
      enddo

      do k = 1,nouter
        zin(k) = xup(k) + eye*yup(k)
      enddo
      call fourierUpsample(nouter,2,zin,zout)
      do k = 1,2*nouter
        xup(k) = dreal(zout(k))
        yup(k) = dimag(zout(k))
      enddo
c     upsample geometry by a factor of 2
      do k = 1,nouter
        zin(k) = denxup(k) + eye*denyup(k)
      enddo
      call fourierUpsample(nouter,2,zin,zout)
      do k = 1,2*nouter
        denxup(k) = dreal(zout(k))
        denyup(k) = dimag(zout(k))
      enddo
c     upsample density by a factor of 2

      call evalLPouter(2*nouter,xup,yup,
     $  denxup,denyup,icount-1,xtar_level,ytar_level,
     $  utar_level,vtar_level,press_tar_level,vort_tar_level)

      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 2 .and. iside(j) .eq. 0) then
          utar(j) = utar(j) + utar_level(icount)
          vtar(j) = vtar(j) + vtar_level(icount)
          press_tar(j) = press_tar(j) + press_tar_level(icount)
          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
          icount = icount + 1
        endif
      enddo
c     END OF DOING POINTS THAT REQUIRE 2X UPSAMPLING


c     START OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 4 .and. iside(j) .eq. 0) then
          xtar_level(icount) = xtar(j)
          ytar_level(icount) = ytar(j)
          icount = icount + 1
        endif
      enddo

      do k = 1,2*nouter
        zin(k) = xup(k) + eye*yup(k)
      enddo
      call fourierUpsample(2*nouter,2,zin,zout)
      do k = 1,4*nouter
        xup(k) = dreal(zout(k))
        yup(k) = dimag(zout(k))
      enddo
c     upsample geometry by a factor of 2
      do k = 1,2*nouter
        zin(k) = denxup(k) + eye*denyup(k)
      enddo
      call fourierUpsample(2*nouter,2,zin,zout)
      do k = 1,4*nouter
        denxup(k) = dreal(zout(k))
        denyup(k) = dimag(zout(k))
      enddo
c     upsample density by a factor of 2

      call evalLPouter(4*nouter,xup,yup,
     $  denxup,denyup,icount-1,xtar_level,ytar_level,
     $  utar_level,vtar_level,press_tar_level,vort_tar_level)

      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 4 .and. iside(j) .eq. 0) then
          utar(j) = utar(j) + utar_level(icount)
          vtar(j) = vtar(j) + vtar_level(icount)
          press_tar(j) = press_tar(j) + press_tar_level(icount)
          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
          icount = icount + 1
        endif
      enddo
c     END OF DOING POINTS THAT REQUIRE 4X UPSAMPLING


c     START OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 8 .and. iside(j) .eq. 0) then
          xtar_level(icount) = xtar(j)
          ytar_level(icount) = ytar(j)
          icount = icount + 1
        endif
      enddo

      do k = 1,4*nouter
        zin(k) = xup(k) + eye*yup(k)
      enddo
      call fourierUpsample(4*nouter,2,zin,zout)
      do k = 1,8*nouter
        xup(k) = dreal(zout(k))
        yup(k) = dimag(zout(k))
      enddo
c     upsample geometry by a factor of 2
      do k = 1,4*nouter
        zin(k) = denxup(k) + eye*denyup(k)
      enddo
      call fourierUpsample(4*nouter,2,zin,zout)
      do k = 1,8*nouter
        denxup(k) = dreal(zout(k))
        denyup(k) = dimag(zout(k))
      enddo
c     upsample density by a factor of 2

      call evalLPouter(8*nouter,xup,yup,
     $  denxup,denyup,icount-1,xtar_level,ytar_level,
     $  utar_level,vtar_level,press_tar_level,vort_tar_level)

      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 8 .and. iside(j) .eq. 0) then
          utar(j) = utar(j) + utar_level(icount)
          vtar(j) = vtar(j) + vtar_level(icount)
          press_tar(j) = press_tar(j) + press_tar_level(icount)
          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
          icount = icount + 1
        endif
      enddo
c     END OF DOING POINTS THAT REQUIRE 8X UPSAMPLING

c     START OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 16 .and. iside(j) .eq. 0) then
          xtar_level(icount) = xtar(j)
          ytar_level(icount) = ytar(j)
          icount = icount + 1
        endif
      enddo

      do k = 1,8*nouter
        zin(k) = xup(k) + eye*yup(k)
      enddo
      call fourierUpsample(8*nouter,2,zin,zout)
      do k = 1,16*nouter
        xup(k) = dreal(zout(k))
        yup(k) = dimag(zout(k))
      enddo
c     upsample geometry by a factor of 2
      do k = 1,8*nouter
        zin(k) = denxup(k) + eye*denyup(k)
      enddo
      call fourierUpsample(8*nouter,2,zin,zout)
      do k = 1,16*nouter
        denxup(k) = dreal(zout(k))
        denyup(k) = dimag(zout(k))
      enddo
c     upsample density by a factor of 2

      call evalLPouter(16*nouter,xup,yup,
     $  denxup,denyup,icount-1,xtar_level,ytar_level,
     $  utar_level,vtar_level,press_tar_level,vort_tar_level)

      icount = 1
      do j=1,ntargets
        if (iup(j) .eq. 16 .and. iside(j) .eq. 0) then
          utar(j) = utar(j) + utar_level(icount)
          vtar(j) = vtar(j) + vtar_level(icount)
          press_tar(j) = press_tar(j) + press_tar_level(icount)
          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
          icount = icount + 1
        endif
      enddo
c     END OF DOING POINTS THAT REQUIRE 16X UPSAMPLING


c     Contribution from Rotlets and Stokeslets
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)
        do j = 1,ntargets
          rx = xtar(j) - centerx(ibod)
          ry = ytar(j) - centery(ibod)
          rho2 = rx**2.d0 + ry**2.d0
          rdots = rx*sto1 + ry*sto2
          utar(j) = utar(j) + 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
          vtar(j) = vtar(j) + 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
          press_tar(j) = press_tar(j) + 1.d0/twopi*
     $        rdots/rho2
          vort_tar(j) = vort_tar(j) + (ry*sto1 - rx*sto2)/
     $        rho2/twopi
c         stokeslet contribution

          utar(j) = utar(j) + rot*ry/rho2
          vtar(j) = vtar(j) - rot*rx/rho2
c         rotlet contribution
        enddo
      enddo
c     Add in contribution from Rotlets and Stokeslets


      do ibod = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (ibod-1)*2*ninner + k)
          deny(k) = den(2*nouter + (ibod-1)*2*ninner + k + ninner)
        enddo

        do j = 1,ninner
          xup(j) = x((ibod-1)*ninner + j)
          yup(j) = y((ibod-1)*ninner + j)
          denxup(j) = denx(j)
          denyup(j) = deny(j)
        enddo
c
        call classifyPoints(ninner,xup,yup,
     $    ntargets,xtar,ytar,iside,nnear,iup)

c       START OF DOING POINTS THAT REQUIRE NO UPSAMPLING
        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 1) then
            xtar_level(icount) = xtar(j)
            ytar_level(icount) = ytar(j)
            icount = icount + 1
          endif
        enddo

        call evalLPbodies(ninner,xup,yup,
     $    denxup,denyup,icount-1,xtar_level,ytar_level,
     $    utar_level,vtar_level,press_tar_level,vort_tar_level)

        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 1) then
            utar(j) = utar(j) + utar_level(icount)
            vtar(j) = vtar(j) + vtar_level(icount)
            press_tar(j) = press_tar(j) + press_tar_level(icount)
            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
            icount = icount + 1
          endif
        enddo
c       END OF DOING POINTS THAT REQUIRE NO UPSAMPLING

c       START OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 2) then
            xtar_level(icount) = xtar(j)
            ytar_level(icount) = ytar(j)
            icount = icount + 1
          endif
        enddo

        do k = 1,ninner
          zin(k) = xup(k) + eye*yup(k)
        enddo
        call fourierUpsample(ninner,2,zin,zout)
        do k = 1,2*ninner
          xup(k) = dreal(zout(k))
          yup(k) = dimag(zout(k))
        enddo
c       upsample geometry by a factor of 2
        do k = 1,ninner
          zin(k) = denxup(k) + eye*denyup(k)
        enddo
        call fourierUpsample(ninner,2,zin,zout)
        do k = 1,2*ninner
          denxup(k) = dreal(zout(k))
          denyup(k) = dimag(zout(k))
        enddo
c       upsample density by a factor of 2

        call evalLPbodies(2*ninner,xup,yup,
     $    denxup,denyup,icount-1,xtar_level,ytar_level,
     $    utar_level,vtar_level,press_tar_level,vort_tar_level)

        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 2) then
            utar(j) = utar(j) + utar_level(icount)
            vtar(j) = vtar(j) + vtar_level(icount)
            press_tar(j) = press_tar(j) + press_tar_level(icount)
            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
            icount = icount + 1
          endif
        enddo
c       END OF DOING POINTS THAT REQUIRE 2X UPSAMPLING

c       START OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 4) then
            xtar_level(icount) = xtar(j)
            ytar_level(icount) = ytar(j)
            icount = icount + 1
          endif
        enddo

        do k = 1,2*ninner
          zin(k) = xup(k) + eye*yup(k)
        enddo
        call fourierUpsample(2*ninner,2,zin,zout)
        do k = 1,4*ninner
          xup(k) = dreal(zout(k))
          yup(k) = dimag(zout(k))
        enddo
c       upsample geometry by a factor of 2
        do k = 1,2*ninner
          zin(k) = denxup(k) + eye*denyup(k)
        enddo
        call fourierUpsample(2*ninner,2,zin,zout)
        do k = 1,4*ninner
          denxup(k) = dreal(zout(k))
          denyup(k) = dimag(zout(k))
        enddo
c       upsample density by a factor of 2

        call evalLPbodies(4*ninner,xup,yup,
     $    denxup,denyup,icount-1,xtar_level,ytar_level,
     $    utar_level,vtar_level,press_tar_level,vort_tar_level)

        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 4) then
            utar(j) = utar(j) + utar_level(icount)
            vtar(j) = vtar(j) + vtar_level(icount)
            press_tar(j) = press_tar(j) + press_tar_level(icount)
            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
            icount = icount + 1
          endif
        enddo
c       END OF DOING POINTS THAT REQUIRE 4X UPSAMPLING

c       START OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 8) then
            xtar_level(icount) = xtar(j)
            ytar_level(icount) = ytar(j)
            icount = icount + 1
          endif
        enddo

        do k = 1,4*ninner
          zin(k) = xup(k) + eye*yup(k)
        enddo
        call fourierUpsample(4*ninner,2,zin,zout)
        do k = 1,8*ninner
          xup(k) = dreal(zout(k))
          yup(k) = dimag(zout(k))
        enddo
c       upsample geometry by a factor of 2
        do k = 1,4*ninner
          zin(k) = denxup(k) + eye*denyup(k)
        enddo
        call fourierUpsample(4*ninner,2,zin,zout)
        do k = 1,8*ninner
          denxup(k) = dreal(zout(k))
          denyup(k) = dimag(zout(k))
        enddo
c       upsample density by a factor of 2

        call evalLPbodies(8*ninner,xup,yup,
     $    denxup,denyup,icount-1,xtar_level,ytar_level,
     $    utar_level,vtar_level,press_tar_level,vort_tar_level)

        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 8) then
            utar(j) = utar(j) + utar_level(icount)
            vtar(j) = vtar(j) + vtar_level(icount)
            press_tar(j) = press_tar(j) + press_tar_level(icount)
            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
            icount = icount + 1
          endif
        enddo
c       END OF DOING POINTS THAT REQUIRE 8X UPSAMPLING

c       START OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 16) then
            xtar_level(icount) = xtar(j)
            ytar_level(icount) = ytar(j)
            icount = icount + 1
          endif
        enddo

        do k = 1,8*ninner
          zin(k) = xup(k) + eye*yup(k)
        enddo
        call fourierUpsample(8*ninner,2,zin,zout)
        do k = 1,16*ninner
          xup(k) = dreal(zout(k))
          yup(k) = dimag(zout(k))
        enddo
c       upsample geometry by a factor of 2
        do k = 1,8*ninner
          zin(k) = denxup(k) + eye*denyup(k)
        enddo
        call fourierUpsample(8*ninner,2,zin,zout)
        do k = 1,16*ninner
          denxup(k) = dreal(zout(k))
          denyup(k) = dimag(zout(k))
        enddo
c       upsample density by a factor of 2

        call evalLPbodies(16*ninner,xup,yup,
     $    denxup,denyup,icount-1,xtar_level,ytar_level,
     $    utar_level,vtar_level,press_tar_level,vort_tar_level)

        icount = 1
        do j=1,ntargets
          if (iup(j) .eq. 16) then
            utar(j) = utar(j) + utar_level(icount)
            vtar(j) = vtar(j) + vtar_level(icount)
            press_tar(j) = press_tar(j) + press_tar_level(icount)
            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
            icount = icount + 1
          endif
        enddo
c       END OF DOING POINTS THAT REQUIRE 16X UPSAMPLING

        do j=1,ntargets
          if (iside(j) .eq. 0 .or. iup(j) .gt. 16) then
            utar(j) = 1.0d20
            vtar(j) = 1.0d20
            press_tar(j) = 1.0d20
            vort_tar(j) = 1.0d20
          endif
        enddo
      enddo

      do k = 1,ntargets
        if (abs(utar(k)) .ge. 1.0d8) then
          utar(k) = 0.d0
          vtar(k) = 0.d0
          press_tar(k) = 0.d0
          vort_tar(k) = 0.d0
        endif
      enddo


      end

c***********************************************************************
c     TODO: THIS CODE SHOULD ONLY EVER BE CALLED IF USING TRAPEZOID RULE
      subroutine evalLPbodies(ninner,x,y,
     $    denx,deny,ntargets,xtar,ytar,
     $    utar,vtar,press_tar,vort_tar)
      implicit real*8 (a-h,o-z)

c     input variables
      dimension x(ninner),y(ninner)
      dimension denx(ninner),deny(ninner)
      dimension xtar(ntargets),ytar(ntargets)

c     output variables
      dimension utar(ntargets),vtar(ntargets)
      dimension press_tar(ntargets),vort_tar(ntargets)

c     local variables
      dimension px(ninner),py(ninner)
      dimension cur(ninner),speed(ninner)
      dimension centerx(1),centery(1)

      call inner_geometry(ninner,1,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do j = 1,ntargets
        utar(j) = 0.d0
        vtar(j) = 0.d0
        press_tar(j) = 0.d0
        vort_tar(j) = 0.d0
        do k = 1,ninner
          rx = xtar(j) - x(k)
          ry = ytar(j) - y(k)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px(k) + ry*py(k)
          rdotden = rx*denx(k) + ry*deny(k)
          dendotn = px(k)*denx(k) +  py(k)*deny(k)
          routn = -rx*py(k) + ry*px(k)
          routden = -rx*deny(k) + ry*denx(k)

          utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $      speed(k)*twopi/dble(ninner)/pi
          vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $      speed(k)*twopi/dble(ninner)/pi
          press_tar(j) = press_tar(j) + 
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed(k)*twopi/dble(ninner)/pi
          vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
     $      rho2**2.d0*speed(k)*twopi/dble(ninner)/pi 
        enddo
      enddo


      end

c***********************************************************************
c     TODO: THIS CODE SHOULD ONLY EVER BE CALLED IF USING TRAPEZOID RULE
      subroutine evalLPouter(n,x,y,
     $    denx,deny,ntargets,xtar,ytar,
     $    utar,vtar,press_tar,vort_tar)
      implicit real*8 (a-h,o-z)

c     input variables
      dimension x(n),y(n)
      dimension denx(n),deny(n)
      dimension xtar(ntargets),ytar(ntargets)

c     output variables
      dimension utar(ntargets),vtar(ntargets)
      dimension press_tar(ntargets),vort_tar(ntargets)

c     local variables
      dimension px(n),py(n)
      dimension cur(n),speed(n)

      call outer_geometry(n,x,y,px,py,cur,speed)
c     build inner geometry

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do j = 1,ntargets
        utar(j) = 0.d0
        vtar(j) = 0.d0
        press_tar(j) = 0.d0
        vort_tar(j) = 0.d0
        do k = 1,n
          rx = xtar(j) - x(k)
          ry = ytar(j) - y(k)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px(k) + ry*py(k)
          rdotden = rx*denx(k) + ry*deny(k)
          dendotn = px(k)*denx(k) +  py(k)*deny(k)
          routn = -rx*py(k) + ry*px(k)
          routden = -rx*deny(k) + ry*denx(k)

          utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $      speed(k)*twopi/dble(n)/pi
          vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $      speed(k)*twopi/dble(n)/pi
          press_tar(j) = press_tar(j) + 
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed(k)*twopi/dble(n)/pi
          vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
     $      rho2**2.d0*speed(k)*twopi/dble(n)/pi 
        enddo
      enddo


      end
c***********************************************************************
      subroutine computeVelocityPressureTargets(ninner,nbodies,nouter,
     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar)
c     Compute the velocity and pressure at a set of target points xtar
c     and ytar.  Target points must be sufficiently far away from the
c     each boundary since no near-singular integration is used.  It is
c     just the vanilla trapezoid rule
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xouter(nouter),youter(nouter)
c     x and y coordinates of confining wall
      dimension px0(nouter), py0(nouter)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nouter), speed0(nouter)
c     Jacobian and curvature of the confining wall
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))

      dimension xtar(ntargets),ytar(ntargets)
      dimension utar(ntargets),vtar(ntargets)
      dimension press_tar(ntargets)

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     build outer geometry

      do j = 1,ntargets
        utar(j) = 0.d0
        vtar(j) = 0.d0
        press_tar(j) = 0.d0
      enddo

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

c     Contribution from the density function on the solid wall
      do k=1,nouter
        denx(k) = den(k)
        deny(k) = den(k+nouter)
      enddo
c     Density function defined on the outer geometry

      do j = 1,ntargets
        do k = 1,nouter
          rx = xtar(j) - xouter(k)
          ry = ytar(j) - youter(k)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px0(k) + ry*py0(k)
          rdotden = rx*denx(k) + ry*deny(k)
          dendotn = px0(k)*denx(k) + py0(k)*deny(k)

          utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $      speed0(k)*twopi/dble(nouter)/pi
          vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $      speed0(k)*twopi/dble(nouter)/pi
          press_tar(j) = press_tar(j) + 
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed0(k)*twopi/dble(nouter)/pi
        enddo
      enddo

c     Contribution from the density function on each inner obstacle
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
        enddo
c       Density function defined on the inner geometry

        do j = 1,ntargets
          do k = 1,ninner
            rx = xtar(j) - x((isou-1)*ninner + k)
            ry = ytar(j) - y((isou-1)*ninner + k)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner + k) + 
     $              ry*py((isou-1)*ninner + k)
            rdotden = rx*denx(k) + ry*deny(k)
            dendotn = px((isou-1)*ninner + k)*denx(k) + 
     $                py((isou-1)*ninner + k)*deny(k)

            utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $        speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
            vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $        speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
            press_tar(j) = press_tar(j) + 
     $          (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $          speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
          enddo
        enddo
      enddo

c     Contribution from Rotlets and Stokeslets
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)
        do j = 1,ntargets
          rx = xtar(j) - centerx(ibod)
          ry = ytar(j) - centery(ibod)
          rho2 = rx**2.d0 + ry**2.d0
          rdots = rx*sto1 + ry*sto2
          utar(j) = utar(j) + 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
          vtar(j) = vtar(j) + 5.d-1/twopi*
     $        (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
          press_tar(j) = press_tar(j) + 1.d0/twopi*
     $        rdots/rho2
c         stokeslet contribution

          utar(j) = utar(j) + rot*ry/rho2
          vtar(j) = vtar(j) - rot*rx/rho2
c         rotlet contribution
        enddo
      enddo
c     Add in contribution from Rotlets and Stokeslets


      end

c***********************************************************************
      subroutine computeVorticityTargets(ninner,nbodies,nouter,
     $    x,y,den,ntargets,xtar,ytar,vort_tar)
c     Compute the vorticity at a set of target points xtar and ytar.  
c     Target points must be sufficiently far away from the each 
c     boundary since no near-singular integration is used.  It is just 
c     the vanilla trapezoid rule
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xouter(nouter),youter(nouter)
c     x and y coordinates of confining wall
      dimension px0(nouter), py0(nouter)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nouter), speed0(nouter)
c     Jacobian and curvature of the confining wall
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))

      dimension xtar(ntargets),ytar(ntargets)
      dimension vort_tar(ntargets)

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     build outer geometry

      do j = 1,ntargets
        vort_tar(j) = 0.d0
      enddo

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

c     Contribution from the density function on the solid wall
      do k=1,nouter
        denx(k) = den(k)
        deny(k) = den(k+nouter)
      enddo
c     Density function defined on the outer geometry

      do j = 1,ntargets
        do k = 1,nouter
          rx = xtar(j) - xouter(k)
          ry = ytar(j) - youter(k)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px0(k) + ry*py0(k)
          rdotden = rx*denx(k) + ry*deny(k)
          routn = -rx*py0(k) + ry*px0(k)
          routden = -rx*deny(k) + ry*denx(k)

          vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
     $      rho2**2.d0*speed0(k)*twopi/dble(nouter)/pi
        enddo
      enddo

c     Contribution from the density function on each inner obstacle
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
        enddo
c       Density function defined on the inner geometry

        do j = 1,ntargets
          do k = 1,ninner
            rx = xtar(j) - x((isou-1)*ninner + k)
            ry = ytar(j) - y((isou-1)*ninner + k)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner + k) + 
     $              ry*py((isou-1)*ninner + k)
            rdotden = rx*denx(k) + ry*deny(k)
            routn = -rx*py((isou-1)*ninner+k) + ry*px((isou-1)*ninner+k)
            routden = -rx*deny(k) + ry*denx(k)

            vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
     $        rho2**2.d0*speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi 
          enddo
        enddo
      enddo

c     Contribution from Rotlets and Stokeslets.  Rotlets have no
c     vorticity
      do ibod = 1,nbodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        do j = 1,ntargets
          rx = xtar(j) - centerx(ibod)
          ry = ytar(j) - centery(ibod)
          rho2 = rx**2.d0 + ry**2.d0
          vort_tar(j) = vort_tar(j) + (ry*sto1 - rx*sto2)/
     $        rho2/twopi
c         stokeslet contribution
        enddo
      enddo
c     Add in contribution from Rotlets and Stokeslets


      end

c***********************************************************************
      subroutine classifyPointsOld(ninner,nbodies,x,y,
     $      ntargets,xtar,ytar,iside,inear)
c     construct a matrix iside which is 1 if a point is inside the fluid
c     domain and 0 if it is outside the fluid domain (ie. in a body).
c     Also construct a matrix inear which is 1 if a point is close to a
c     boundary and 0 otherwise

      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xtar(ntargets),ytar(ntargets)
      dimension iside(ntargets),inear(ntargets)
      
      dimension arclengths(nbodies)
      dimension dist_bodies(nbodies)

      twopi = 8.d0*datan(1.d0)

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      do j = 1,nbodies
        arclengths(j) = 0.d0
        do k = (j-1)*ninner+1,j*ninner
          arclengths(j) = arclengths(j) + speed(k)
        enddo
        arclengths(j) = arclengths(j)*twopi/dble(ninner)
        arclengths(j) = arclengths(j)/dble(ninner)
c       ds term
      enddo

      do i = 1,ntargets
        inear(i) = 0.d0
c       assume point is not too close to a boundary
        ibodycp = 1
        indexcp = 1
c       index of the closest body to target point and the index
c       of the closest point on that body

        dist_min = 1.d10
        do j = 1,nbodies
          dist_bodies(j) = 1.d10
          do k = 1,ninner
            ind = (j-1)*ninner + k
            dist2 = (xtar(i) - x(ind))**2.d0 + (ytar(i) - y(ind))**2.d0
            if (dist2 .lt. dist_bodies(j)) then
              dist_bodies(j) = dist2
            endif
            if (dist2 .lt. dist_min) then
              dist_min = dist2
              ibodycp = j
              indexcp = k
            endif
          enddo

          if (dsqrt(dist_bodies(j)) .lt. 5.d0*arclengths(j)) then
            inear(i) = 1
          endif
        enddo

        ind = (ibodycp - 1)*ninner + indexcp
        vec_dot_n = (xtar(i) - x(ind))*px(ind) + 
     $              (ytar(i) - y(ind))*py(ind)
c       This quanitity will be negative if it is a point in the fluid
c       and it will be positive if it is a point in a body

        if (vec_dot_n .lt. 0.d0) then
          iside(i) = 1
        else
          iside(i) = 0
        endif


      enddo
            

      end

c***********************************************************************
      subroutine classifyPoints(ninner,x,y,
     $      ntargets,xtar,ytar,iside,nnear,iup)
c     construct a matrix iside which is 1 if a point is inside the fluid
c     domain and 0 if it is outside the fluid domain (ie. in a body).
c     Also construct a matrix inear which is 1 if a point is close to a
c     boundary and 0 otherwise

      implicit real*8 (a-h,o-z)

      dimension x(ninner),y(ninner)
      dimension px(ninner),py(ninner)
      dimension cur(ninner),speed(ninner)
      dimension centerx(1),centery(1)
      dimension xtar(ntargets),ytar(ntargets)
      dimension iside(ntargets)
      dimension iup(ntargets)

      complex *16 z(ninner),dz(ninner),d2z(ninner)
      complex *16 h,dh,d2h
      dimension wsave(4*ninner + 15)
      complex *16 eye
      dimension modes(ninner)
      
      eye = (0.d0,1.d0)
      twopi = 8.d0*datan(1.d0)

      do k = 1,ninner/2
        modes(k) = k-1
      enddo
      do k=ninner/2+1,ninner
        modes(k) = -ninner + k - 1
      enddo

      call inner_geometry(ninner,1,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      call DCFFTI(ninner,wsave)
      do k = 1,ninner
        z(k) = x(k) + eye*y(k)
      enddo

      call DCFFTF(ninner,z,wsave)
      do k = 1,ninner
        dz(k) = eye*modes(k)*z(k)
        d2z(k) = -modes(k)*modes(k)*z(k)
      enddo
c     Compute shape and its first two derivatives in Fourier space

      nnear = 0
      arclength = 0.d0
      do k = 1,ninner
        arclength = arclength + speed(k)
      enddo
      arclength = arclength*twopi/dble(ninner)
      arclength = arclength/dble(ninner)
c     ds term

      do i = 1,ntargets
        indexcp = 1
c       index of the closest body to target point and the index
c       of the closest point on that body

        dist_body = 1.d10
        do k = 1,ninner
          dist2 = (xtar(i) - x(k))**2.d0 + (ytar(i) - y(k))**2.d0
          if (dist2 .lt. dist_body) then
            dist_body = dist2
            indexcp = k
          endif
        enddo
c       closest discretization point, but this isn't good enough to 
c       use to define what upsampling rate is necessary.  Need to do
c       a few Newton iterations to find the actual distance.

        theta = dble(indexcp-1)*twopi/dble(ninner)
        do inewton = 1,5
          h = 0.d0
          dh = 0.d0
          d2h = 0.d0
          do k = 1,ninner
            h = h + z(k)*exp(eye*modes(k)*theta)
            dh = dh + dz(k)*exp(eye*modes(k)*theta)
            d2h = d2h + d2z(k)*exp(eye*modes(k)*theta)
          enddo
          h = h/dble(ninner)
          dh = dh/dble(ninner)
          d2h = d2h/dble(ninner)
          if (dabs(centerx(1) - 0.5d0) .le. 1.d-8 
     $          .and. i .eq. 7734) then
          endif

          f = (dreal(h) - xtar(i))*dreal(dh) + 
     $        (dimag(h) - ytar(i))*dimag(dh)
          df = (dreal(h) - xtar(i))*dreal(d2h) + 
     $         (dimag(h) - ytar(i))*dimag(d2h) + 
     $         abs(dh)**2.d0

          theta = theta - f/df
        enddo

        h = 0.d0
        do k = 1,ninner
          h = h + z(k)*exp(eye*modes(k)*theta)
        enddo
        h = h/dble(ninner)
        dist_body = (dreal(h) - xtar(i))**2.d0 + 
     $              (dimag(h) - ytar(i))**2.d0

        if (dsqrt(dist_body) .lt. 3.125d-1*arclength) then
          iup(i) = 100
        elseif (dsqrt(dist_body) .lt. 6.25d-1*arclength) then
          iup(i) = 16
        elseif (dsqrt(dist_body) .lt. 1.25d0*arclength) then
          iup(i) = 8 
        elseif (dsqrt(dist_body) .lt. 2.5d0*arclength) then
          iup(i) = 4 
        elseif (dsqrt(dist_body) .lt. 5.0d0*arclength) then
          iup(i) = 2 
        elseif (dsqrt(dist_body) .ge. 5.0d0*arclength) then
          iup(i) = 1 
        endif


        vec_dot_n = (xtar(i) - x(indexcp))*px(indexcp) + 
     $              (ytar(i) - y(indexcp))*py(indexcp)
c       This quanitity will be negative if it is a point in the fluid
c       and it will be positive if it is a point in a body

        if (vec_dot_n .lt. 0.d0) then
          iside(i) = 1
        else
          iside(i) = 0
        endif

      enddo


      end

c***********************************************************************
      subroutine computePressure(ninner,nbodies,nouter,x,y,den,
     $    pressure)
c     Compute the pressure on the boundary of each grain.  Need to 
c     remove one order of the singularity by multiplying by the
c     density function at the target point and then use odd-even 
c     integration
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
c     x and y coordinates of the normal of the obstacle
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
c     Jacobian and curvature of the geometry
      dimension centerx(nbodies),centery(nbodies)

      dimension xouter(nouter),youter(nouter)
c     x and y coordinates of confining wall
      dimension px0(nouter), py0(nouter)
c     x and y coordinates of the normal of the confining wall
      dimension cur0(nouter), speed0(nouter)
c     Jacobian and curvature of the confining wall

      dimension den(2*ninner*nbodies + 3*nbodies + 2*nouter)
c     x and y coordinates of the density function

      dimension pressure(ninner*nbodies)
c     shear stress on the obstacle

      dimension denx(max(ninner,nouter))
      dimension deny(max(ninner,nouter))
      complex *16 zden(ninner)
      dimension wsave(4*ninner + 15)
      complex *16 eye

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     build outer geometry

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      eye = (0.d0,1.d0)

      do j = 1,ninner*nbodies
        pressure(j) = 0.d0
      enddo

c     START OF SOURCE POINTS == SOLID WALL, TARGET POINTS == OBSTACLES
      do k=1,nouter
        denx(k) = den(k)
        deny(k) = den(k+nouter)
      enddo
c     Density function defined on the outer geometry

      do itar = 1,nbodies
c       loop over target bodies
        do k = 1,ninner
c         loop over target points
          do j = 1,nouter
c           loop over sources
            rx = x((itar-1)*ninner + k) - xouter(j)
            ry = y((itar-1)*ninner + k) - youter(j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px0(j) + ry*py0(j)
            rdotden = rx*denx(j) + ry*deny(j)
            dendotn = px0(j)*denx(j) + py0(j)*deny(j)

            pressure((itar-1)*ninner + k) = 
     $        pressure((itar-1)*ninner + k) + 
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed0(j)*twopi/dble(nouter)/pi
          enddo
        enddo
      enddo
c     END OF SOURCE POINTS == SOLID WALL, TARGET POINTS == OBSTACLES

c     START OF SOURCE POINTS == OBSTACLES, 
c              TARGET POINTS == OBSTACLES
      do isou = 1,nbodies
        do j = 1,ninner
          denx(j) = den(2*nouter + (isou-1)*2*ninner + j)
          deny(j) = den(2*nouter + (isou-1)*2*ninner + j + ninner)
        enddo
c       density function due to obstacle isou

c       START OF TARGET POINTS == OBSTACLE isou
c       loop over odd-inexed target points
        do k = 1,ninner,2
c         loop over source points
          do j = 2,ninner,2
            sx = denx(j) - denx(k)
            sy = deny(j) - deny(k)
c           subtract off the density at the target point to reduce
c           the strenght of the singularity to 1/r.  Then, odd-even
c           integration can be applied and is guaranteed to converge

            rx = x((isou-1)*ninner+k) - x((isou-1)*ninner+j)
            ry = y((isou-1)*ninner+k) - y((isou-1)*ninner+j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner+j) +
     $              ry*py((isou-1)*ninner+j)
            rdotden = rx*sx + ry*sy
            dendotn = px((isou-1)*ninner+j)*sx + 
     $                py((isou-1)*ninner+j)*sy
            
            pressure((isou-1)*ninner+k) =
     $        pressure((isou-1)*ninner+k) + 2.d0*
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed((isou-1)*ninner + j)*twopi/dble(ninner)/pi
c           need to multiply by 2 since the grid spacing is twice as
c           large because of the odd-even integration
          enddo
        enddo

c       loop over even-indexed target points
        do k = 2,ninner,2
c         loop over source points
          do j = 1,ninner,2
            sx = denx(j) - denx(k)
            sy = deny(j) - deny(k)
c           subtract off the density at the target point to reduce
c           the strenght of the singularity to 1/r.  Then, odd-even
c           integration can be applied and is guaranteed to converge

            rx = x((isou-1)*ninner+k) - x((isou-1)*ninner+j)
            ry = y((isou-1)*ninner+k) - y((isou-1)*ninner+j)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px((isou-1)*ninner+j) +
     $              ry*py((isou-1)*ninner+j)
            rdotden = rx*sx + ry*sy
            dendotn = px((isou-1)*ninner+j)*sx + 
     $                py((isou-1)*ninner+j)*sy
            
            pressure((isou-1)*ninner+k) =
     $        pressure((isou-1)*ninner+k) + 2.d0*
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed((isou-1)*ninner + j)*twopi/dble(ninner)/pi
c           need to multiply by 2 since the grid spacing is twice as
c           large because of the odd-even integration
          enddo
        enddo
c       END OF TARGET POINTS == OBSTACLE isou

c       START OF TARGET POINTS ~= OBSTACLE isou
        do itar = 1,nbodies
          if (itar .eq. isou) then
            cycle
          endif
c         skip the diagonal term since this was taking care of above
c         with the trapezoid rule with the correcting limiting value at 
c         the diagonal

c         loop over target points
          do k = 1,ninner

c           loop over source points
            do j = 1,ninner
              rx = x((itar-1)*ninner+k) - x((isou-1)*ninner+j)
              ry = y((itar-1)*ninner+k) - y((isou-1)*ninner+j)
              rho2 = rx**2.d0 + ry**2.d0
              rdotn = rx*px((isou-1)*ninner+j) +
     $                ry*py((isou-1)*ninner+j)
              rdotden = rx*denx(j) + ry*deny(j)
              dendotn = px((isou-1)*ninner+j)*denx(j) + 
     $                  py((isou-1)*ninner+j)*deny(j)

              pressure((itar-1)*ninner+k) =
     $          pressure((itar-1)*ninner+k) + 
     $          (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)* 
     $          speed((isou-1)*ninner + j)*twopi/dble(ninner)/pi
            enddo
          enddo
        enddo
c       END OF TARGET POINTS ~= OBSTACLE isou

      enddo
c     END OF SOURCE POINTS == OBSTACLES, 
c            TARGET POINTS == OBSTACLES


c     START OF JUMP ALONG DIAGONAL TERM
      call DCFFTI(ninner,wsave)
      do ibod = 1,nbodies
        do k = 1,ninner
          zden(k) = den(2*nouter+(ibod-1)*2*ninner+k) +
     $          eye*den(2*nouter+(ibod-1)*2*ninner+k+ninner)
        enddo
        call fourierDiff(ninner,zden,wsave)
c       real part of zden is parameter derivative of the first 
c       component of the density function, and the imaginary part is
c       the derivative of the second component of the density function

        do k = 1,ninner
          tx = -py((ibod-1)*ninner+k)
          ty =  px((ibod-1)*ninner+k)
c         tangent vector
          dsdtx = dreal(zden(k))/speed((ibod-1)*ninner+k)
c         derivative of the first component with respect to arclength
          dsdty = dimag(zden(k))/speed((ibod-1)*ninner+k)
c         derivative of the second component with respect to arclength
          dsdt_dot_tau = dsdtx*tx + dsdty*ty
c         dot product of the tangential derivative of the density
c         function with the tangent vector
          pressure((ibod-1)*ninner+k) = pressure((ibod-1)*ninner+k) -
     $        1.d0*dsdt_dot_tau
c         add in jump condition
        enddo
      enddo
c     END OF JUMP ALONG DIAGONAL TERM


c     START OF SOURCE POINTS == ROTLETS AND STOKESLETS, 
c              TARGET POINTS == OBSTACLES
      do ibod = 1,nbodies
c       loop over target bodies
        sto1 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+1)
        sto2 = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+2)
        rot  = den(2*nouter+2*ninner*nbodies+(ibod-1)*3+3)
c       Stokeslet flow is not pressure-free
c       Rotlet flow is pressure-free
        do itar = 1,nbodies
c         loop over target points
          do k = 1,ninner
            rx = x((itar-1)*ninner+k) - centerx(ibod)
            ry = y((itar-1)*ninner+k) - centery(ibod)
            rho2 = rx**2.d0 + ry**2.d0
            rdots = rx*sto1 + ry*sto2
            pressure((itar-1)*ninner + k) = 
     $          pressure((itar-1)*ninner + k) + 1.d0/twopi*
     $          rdots/rho2
          enddo
        enddo
      enddo
c     END OF SOURCE POINTS == ROTLETS AND STOKESLETS, 
c            TARGET POINTS == OBSTACLES

      end



c***********************************************************************
      subroutine computeDrag(ninner,nbodies,x,y,
     $    shear_stress,pressure,drag)
c     Compute the drag of each grain
c     THIS IS NOT BEING USED AS THERE IS A JULIA VERSION SO THAT WE WERE
c     CLEAR WITH THE NORMAL AND TANGENTIAL VECTOR CONVENTIONS
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
c     x and y coordinates of the normal of the obstacle
      dimension cur(ninner*nbodies),speed(ninner*nbodies)
c     Jacobian and curvature of the geometry

      dimension shear_stress(ninner*nbodies)
      dimension pressure(ninner*nbodies)
c     shear stress and pressure on the obstacle

      dimension drag(2*nbodies)

      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      twopi = 8.d0*datan(1.d0)

      do k = 1,2*nbodies
        drag(k) = 0.d0
      enddo
c     initialize drag to be zero

      do k = 1,nbodies
        do j = 1,ninner
          drag((k-1)*2 + 1) = drag((k-1)*2 + 1) +
     $      (pressure((k-1)*ninner + j)*px((k-1)*ninner + j) - 
     $      shear_stress((k-1)*ninner + j)*py((k-1)*ninner + j))*
     $      speed((k-1)*ninner + j)
          drag((k-1)*2 + 2) = drag((k-1)*2 + 2) +
     $      (pressure((k-1)*ninner + j)*py((k-1)*ninner + j) + 
     $      shear_stress((k-1)*ninner + j)*px((k-1)*ninner + j))*
     $      speed((k-1)*ninner + j)
        enddo
        drag((k-1)*2 + 1) = drag((k-1)*2 + 1)/twopi
        drag((k-1)*2 + 2) = drag((k-1)*2 + 2)/twopi
      enddo


      end
