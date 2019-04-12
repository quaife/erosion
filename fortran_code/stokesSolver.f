      subroutine stokesSolver(nninner,nnbodies,nnouter,
     $     ifmm,ibary,maxl,xx,yy,den,iter)
c     Input x and y coordinates and return the density function on the
c     boundary.  Outer wall is used to enclose the inner boundary so
c     that Stokes paradox is avoided
      implicit real*8 (a-h,o-z)

      dimension xx(nninner*nnbodies),yy(nninner*nnbodies)
c     x and y coordinates of obstacle

      parameter (nmax = 2**17)
      parameter (maxbodies = 500)
      parameter (ntargets = 2500)
c     max points on the boundary of the obstacle
      parameter (liwork = 30)
!      parameter (maxl = 2000, liwork = 30)
!      parameter (lrwork = 10 + nmax*(maxl+6) + 
!     $     maxl*(maxl+3))
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

      real *8, allocatable :: gmwork(:)
      dimension igwork(liwork)
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
c     cant declare input variables into a common field, so need new
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
      print *, 'load Bd condition'
      
      lrwork = 10 + nmax*(maxl+6) + maxl*(maxl+3)
      allocate(gmwork(lrwork))
      
      call solveBIE(ninner,nbodies,nouter,den,rhs,
     $      gmwork,lrwork,igwork,liwork,maxl,ifmm,ibary,iter)
c     solve for the density function with GMRES
      print *, 'Solve BIE'
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
      subroutine fourierDiff(n,den)!,wsave)
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
        call fourierDiff(n,zden1)!,wsave)
        zden2 = zden1
        call fourierDiff(n,zden2)!,wsave)
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

      external msolve_DLP,matvec_DLP,matvec_DLP_fmm,matvec_DLP_bary,
     $         matvec_DLP_fmmbary 

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

      print *, 'Start GMRES'
      call cpu_time(t1)
      if (ibary .eq. 1) then
        if (ifmm .eq. 1) then
          print*,'USING FMMBarycentric'         
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
     $      nelt,ia,ja,a,isym,
     $      matvec_DLP_fmmbary,msolve_DLP,itol,tol,itmax,iter,err,
     $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
        else
          print*,'USING Barycentric'
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
     $      nelt,ia,ja,a,isym,
     $      matvec_DLP_bary,msolve_DLP,itol,tol,itmax,iter,err,
     $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c       Use Barycentric
        endif
      else
        if (ifmm .eq. 1) then
          print*,'USING FMM'
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
     $      nelt,ia,ja,a,isym,
     $      matvec_DLP_fmm,msolve_DLP,itol,tol,itmax,iter,err,
     $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c        USE FMM
        else
          print*,'USING DIRECT'
          call DGMRES(2*nouter+2*ninner*nbodies+3*nbodies,rhs,den,
     $      nelt,ia,ja,a,isym,
     $      matvec_DLP,msolve_DLP,itol,tol,itmax,iter,err,
     $      ierr,6,sb,sx,gmwork,lrwork,igwork,liwork,rwork,iwork)
c        DON'T USE FMM
        endif
      endif
c     use GMRES to find the density function of the double layer
c     potential
      call cpu_time(t2)
      print *, 'GMRES  time =', t2-t1



      end

c***********************************************************************

      subroutine matvec_DLP_fmmbary(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 500)

      dimension den(ntotal)
      dimension vel(ntotal)

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter

      dimension x(nmax),y(nmax)
      dimension xtar(nbodies*ninner),ytar(nbodies*ninner),
     &          xsou(max(ninner,nouter)),ysou(max(ninner,nouter))
      dimension xtarloc(ninner),ytarloc(ninner)
      dimension px(nmax),py(nmax)
      dimension cur(nmax),speed(nmax)
      dimension centerx(maxbodies),centery(maxbodies)

      dimension xouter(nmax),youter(nmax)
      dimension px0(nmax),py0(nmax)
      dimension cur0(nmax),speed0(nmax)

      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))
      dimension denxb(max(ninner,nouter)),denyb(max(ninner,nouter))       
      dimension ux(max(ninner*nbodies,nouter)),
     &          uy(max(ninner*nbodies,nouter))
      dimension u1x(max(ninner*nbodies,nouter)),
     &          u1y(max(ninner*nbodies,nouter)),
     &          u2x(max(ninner*nbodies,nouter)),
     &          u2y(max(ninner*nbodies,nouter))       
      
      dimension indw2b(2,ninner*nbodies)
      dimension indb2b(2,ninner*nbodies)
      dimension indb2w(nouter)

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
      nbeta = 1

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

c    START OF REPLACING THE TRAPEZOID RUSELTS BY BARYCENTRIC ONE WHEN
c    TWO POINTS ARE TOO CLOSED


c     START the wall-bodies section
c     start the wall-to-bodies part
      do k = 1,nouter
        denx(k) = den(k)
        deny(k) = den(k + nouter)
      enddo      

      m=1
      do itar = 1,nbodies

c       loop over target points
c       when the distancce between wall and bodies is closed
        if( abs(centery(itar)) .ge. 0.5d0) then

          do k=1,ninner
            xtarloc(k)=x((itar-1)*ninner+k)
            ytarloc(k)=y((itar-1)*ninner+k)
            denxb(k) = den(2*nouter + (itar-1)*2*ninner + k)
            denyb(k) = den(2*nouter + (itar-1)*2*ninner + k + ninner)
c           save xtarloc, ytarloc, denxb, and denyb for body-to-wall barycentric           
          
            if(abs(y((itar-1)*ninner+k)) .ge. 0.5d0)  then
              ux(k) = 0.d0
              uy(k) = 0.d0
              xtar(m) = x((itar-1)*ninner+k)
              ytar(m) = y((itar-1)*ninner+k)
              indw2b(1,m) = itar
              indw2b(2,m) = k
              m = m + 1            
c             save xtar, ytar, indw2b for wall-to-bodies barycentric            
              
c             Apply trapezoid rule from wall to point_indw2b
c             loop over source points
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
              
              vel(2*nouter+(itar-1)*2*ninner+k) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k) - ux(k)
              vel(2*nouter+(itar-1)*2*ninner+k+ninner) = 
     $          vel(2*nouter+(itar-1)*2*ninner+k+ninner) - uy(k)
c             Take out the trapezoid part at p_indw2b     
            endif
          enddo
c         end of the wall-to-bodies part and save the index indw2b for doing barycentric 
c         later once for all 

c         start the body-to-wall part
          rtar = 0.d0
          do k = 1, ninner
            rtar = rtar + speed((itar-1)*ninner+k)  
          enddo
          rtar = rtar/ninner
          i = 1
          do j=1,nouter
            dx = centerx(itar) - xouter(j)
            dy = centery(itar) - youter(j)
            d2 = dx**2.d0 + dy**2.d0
            r = 1.d0 - abs(centery(itar)) + 2.d0*rtar
            
            if( d2 .le. r**2) then
              xsou(i) = xouter(j)
              ysou(i) = youter(j)
              indb2w(i) = j              
              i = i + 1
            endif
          enddo
          
          do k=1,i-1
            ux(k) = 0.d0
            uy(k) = 0.d0
          
            do j=1,ninner
              rx = xsou(k) - x((itar-1)*ninner + j)
              ry = ysou(k) - y((itar-1)*ninner + j)
              rdotn = rx*px((itar-1)*ninner + j) + 
     $              ry*py((itar-1)*ninner + j)
              rdotden = rx*denxb(j) + ry*denyb(j)
              rho2 = rx**2.d0 + ry**2.d0

              ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*
     $            speed((itar-1)*ninner+j)
              uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*
     $            speed((itar-1)*ninner+j)
            enddo
            ux(k) = ux(k)*twopi/dble(ninner)/pi
            uy(k) = uy(k)*twopi/dble(ninner)/pi
          enddo

          do k=1,i-1
            vel(indb2w(k)) = vel(indb2w(k)) - ux(k)
            vel(indb2w(k) + nouter) = vel(indb2w(k) + nouter) - uy(k)
          enddo
          
          ninnc = ninner/nbeta
          nder = 1
          call StokesExteriorDLP(ninner,ninnc,nbeta,xtarloc,ytarloc,
     &          denxb,denyb,px,py,i-1,xsou,ysou,itar,nbodies,
     &          nder,ux,uy,u1x,u1y,u2x,u2y)
     
          do k=1,i-1
            vel(indb2w(k)) = vel(indb2w(k)) + ux(k)
            vel(indb2w(k) + nouter) = vel(indb2w(k) + nouter) + uy(k)
          enddo           
        endif  
      enddo
c     end of the body-to-wall part

c     Apply barycentric rule  from wall to point_indw2b once for all      
      noutc = nouter/nbeta
      nder = 1
      call StokesInteriorDLP(nouter,noutc,nbeta,xouter,youter,
     &      denx,deny,px0,py0,m-1,xtar,ytar,nder,
     &      ux,uy,u1x,u1y,u2x,u2y)
      
      do k = 1,m-1
        vel(2*nouter+(indw2b(1,k)-1)*2*ninner+indw2b(2,k)) = 
     $    vel(2*nouter+(indw2b(1,k)-1)*2*ninner+indw2b(2,k)) + ux(k)
        vel(2*nouter+(indw2b(1,k)-1)*2*ninner+indw2b(2,k)+ninner) = 
     $    vel(2*nouter+(indw2b(1,k)-1)*2*ninner+indw2b(2,k)+ninner) + 
     $    uy(k)
      enddo      
c     END of the section of the wall-bodies replacement

c     START OF the section of body-to-bodies replacement
      do isou = 1,nbodies
        rsou = 0.d0
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
          xsou(k) = x((isou-1)*ninner+k)
          ysou(k) = y((isou-1)*ninner+k)
          rsou = rsou + speed((isou-1)*ninner+k)
        enddo
        rsou = rsou/ninner
c       density function due to obstacle isou

c       START OF TARGET POINTS ~= OBSTACLE isou
        l = 1
c        print *, 'isou', isou
        do itar = 1,nbodies
          if (itar .eq. isou) then
            cycle
          endif
          rtar = 0.d0
          do k=1,ninner
            rtar = rtar +  speed((itar-1)*ninner+k)
          enddo
          rtar = rtar/ninner
c         skip the diagonal term since this was taking care above with
c         the trapezoid rule with the correcting liming value at the
c         diagonal
          dx = centerx(isou) - centerx(itar)
          dy = centery(isou) - centery(itar)
          d2 = dx**2.d0+dy**2.d0
          rsum = 2.d0*(rsou + rtar)
c         check if B_itar is closed to B_isou          
          if( d2 .le. rsum**2) then 
c         loop over target points
            
            do k = 1,ninner
              ux(k) = 0.d0
              uy(k) = 0.d0
              dx = centerx(isou)-x((itar-1)*ninner+k)
              dy = centery(isou)-y((itar-1)*ninner+k)              
              d2 = dx**2.d0 + dy**2.d0
c             Check if point_k is closed to B_isou              
              if (d2 .le. rsum**2) then              
                xtar(l) = x((itar-1)*ninner+k)
                ytar(l) = y((itar-1)*ninner+k)
                indb2b(1,l) = itar
                indb2b(2,l) = k               
c               Save the index indb2b for appling barycentric from B_sou to 
c               point_indb2b of B_itar if they are closed                
                l = l + 1

                do j=1,ninner
                  rx = x((itar-1)*ninner+k) - x((isou-1)*ninner+j)
                  ry = y((itar-1)*ninner+k) - y((isou-1)*ninner+j)
                  rdotn = rx*px((isou-1)*ninner+j) + 
     $                ry*py((isou-1)*ninner+j) 
                  rdotden = rx*denx(j) + ry*deny(j)
                  rho2 = rx**2.d0 + ry**2.d0

                  ux(k) = ux(k) + rdotn*rdotden/rho2/rho2*rx*
     $              speed((isou-1)*ninner+j)
                  uy(k) = uy(k) + rdotn*rdotden/rho2/rho2*ry*
     $              speed((isou-1)*ninner+j)
                enddo
              ux(k) = ux(k)*twopi/dble(ninner)/pi
              uy(k) = uy(k)*twopi/dble(ninner)/pi
c             Get the result by trapezoid rule at p_indb2b
              vel(2*nouter + (itar-1)*2*ninner + k) = 
     $          vel(2*nouter + (itar-1)*2*ninner + k) - ux(k)
              vel(2*nouter + (itar-1)*2*ninner + k + ninner) = 
     $          vel(2*nouter + (itar-1)*2*ninner + k + ninner) - uy(k)
c             Take the result away at point_indb2b and replace it by 
c             the result of barycentric later              
              endif
            enddo
          endif
        enddo
        if( l .gt. 1) then
c         Apply barycentric rule from B_isou and point_indb2b          
          ninnc = ninner/nbeta
          nder = 1
          call StokesExteriorDLP(ninner,ninnc,nbeta,xsou,ysou,
     &      denx,deny,px,py,l-1,xtar,ytar,isou,nbodies,nder,
     &      ux,uy,u1x,u1y,u2x,u2y)
     
          do k = 1,l-1          
            vel(2*nouter+(indb2b(1,k)-1)*2*ninner+indb2b(2,k)) = 
     $        vel(2*nouter+(indb2b(1,k)-1)*2*ninner+indb2b(2,k)) + 
     $        ux(k)
            vel(2*nouter+(indb2b(1,k)-1)*2*ninner+indb2b(2,k)+ninner) = 
     $      vel(2*nouter+(indb2b(1,k)-1)*2*ninner+indb2b(2,k)+ninner) + 
     $        uy(k)        
          enddo
        endif  
      enddo          

c    END of the section of body to bodies replacement
  

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
      subroutine matvec_DLP_bary(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 500)

      dimension den(ntotal)
      dimension vel(ntotal)

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter

      dimension x(nmax),y(nmax)
      dimension xtar(ninner),ytar(ninner),
     &          xsou(ninner),ysou(ninner) 
      dimension px(nmax),py(nmax)
      dimension cur(nmax),speed(nmax)
      dimension centerx(maxbodies),centery(maxbodies)

      dimension xouter(nmax),youter(nmax)
      dimension px0(nmax),py0(nmax)
      dimension cur0(nmax),speed0(nmax)

      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))
      dimension ux(max(ninner,nouter)),uy(max(ninner,nouter))
      dimension u1x(max(ninner,nouter)),u1y(max(ninner,nouter)),
     &          u2x(max(ninner,nouter)),u2y(max(ninner,nouter))  
      
      dimension uxtar(ninner),uytar(ninner),
     &          ux1tar(ninner),ux2tar(ninner),
     &          uy1tar(ninner),uy2tar(ninner)
      dimension press_tar(ninner)

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
          xtar(k) = x((itar-1)*ninner+k)
          ytar(k) = y((itar-1)*ninner+k)           
        enddo
        
c      loop over source points
       nbeta = 1
       noutc = nouter
       nder = 1
       call StokesInteriorDLP(nouter,noutc,nbeta,xouter,youter,
     &          denx,deny,px0,py0,ninner,xtar,ytar,nder,
     &          ux,uy,u1x,u1y,u2x,u2y) 

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
          xsou(k) = x((isou-1)*ninner+k)
          ysou(k) = y((isou-1)*ninner+k)           
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
            xtar(k) = x((itar-1)*ninner+k)
            ytar(k) = y((itar-1)*ninner+k)           
          enddo
            
c         loop over source points
          nbeta = 1
          ninnc = ninner
          nder = 1
          call StokesExteriorDLP(ninner,ninnc,nbeta,xsou,ysou,
     &          denx,deny,px,py,ninner,xtar,ytar,isou,nbodies,nder,
     &          ux,uy,u1x,u1y,u2x,u2y) 
c        print *, '2'       
c          do k = 1,ninner
c            ux(k) = uxtar(k)
c            uy(k) = uytar(k)
c          enddo
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
        enddo
c       loop over source points
        nbeta = 1
        ninnc = ninner
        nder = 1
        call StokesExteriorDLP(ninner,ninnc,nbeta,xsou,ysou,
     &          denx,deny,px,py,nouter,xouter,youter,isou,nbodies,nder,
     &          ux,uy,u1x,u1y,u2x,u2y) 
c        print *, '3'       
c        do k=1,nouter
c          ux(k) = uxtar(k)
c          uy(k) = uytar(k)
c        enddo
        
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
c        print *, '4'  
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
      subroutine matvec_DLP(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 500)

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


c      do k = 1,2*nouter + 2*ninner*nbodies
c        print*,vel(k)
c      enddo

      return
      end      
c***********************************************************************
      subroutine matvec_DLP_fmm(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**17)
      parameter (maxbodies = 500)

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
      subroutine deformation_on_boundary_trap(ninner,nbodies,x,y,
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
        call fourierDiff(ninner,zden)!,wsave)
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
      subroutine deformation_on_boundary(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor at points on the boundary of the
c     obstacle
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension xtar(ninner),ytar(ninner),
     &          xsou(ninner),ysou(ninner)       
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
      
      real *8 ux(ninner),uy(ninner)
      real *8 u1x(ninner),u1y(ninner),
     &        u2x(ninner),u2y(ninner)        

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

c      trapezoid case
c      do itar=1,nbodies
cc       loop over targets
c        do k=1,ninner
cc         loop over sources
c          do j=1,nouter
c            rx = x((itar-1)*ninner+k) - xouter(j)
c            ry = y((itar-1)*ninner+k) - youter(j)
c            rho2 = rx**2.d0 + ry**2.d0
c            rdotn = rx*px0(j) + ry*py0(j)
c            routn = rx*py0(j) + ry*px0(j)
c            rdotden = rx*denx(j) + ry*deny(j)
c            routden = rx*deny(j) + ry*denx(j)
c
c            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + 5.d-1*
c     $          (2.d0*rdotn*rdotden/rho2/rho2 + 
c     $          2.d0*rdotden/rho2/rho2*rx*px0(j) +
c     $          2.d0*rdotn/rho2/rho2*rx*denx(j) -
c     $          8.d0*rdotn*rdotden/rho2**3.d0*rx*rx)*
c     $          speed0(j)*twopi/dble(nouter)/pi
c
c            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 5.d-1*
c     $          (rdotden*routn/rho2/rho2 +
c     $          rdotn*routden/rho2/rho2 -
c     $          8.d0*rdotn*rdotden*rx*ry/rho2**3.d0)*
c     $          speed0(j)*twopi/dble(nouter)/pi
c
c            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 5.d-1*
c     $          (2.d0*rdotn*rdotden/rho2/rho2 + 
c     $          2.d0*rdotden/rho2/rho2*ry*py0(j) +
c     $          2.d0*rdotn/rho2/rho2*ry*deny(j) -
c     $          8.d0*rdotn*rdotden/rho2**3.d0*ry*ry)*
c     $          speed0(j)*twopi/dble(nouter)/pi
c          enddo
c        enddo
c      enddo

c     barycentric case
      do itar=1,nbodies
c       loop over targets
        do k=1,ninner
          xtar(k) = x((itar-1)*ninner+k)
          ytar(k) = y((itar-1)*ninner+k)           
        enddo
       nder = 2 
       nbeta = 1
       noutc = nouter
       call StokesInteriorDLP(nouter,noutc,nbeta,xouter,youter,
     &          denx,deny,px0,py0,ninner,xtar,ytar,nder,
     &          ux,uy,u1x,u1y,u2x,u2y)
       
        
        do k = 1,ninner
            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + u1x(k) 
            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 
     &                               0.5d0*(u1y(k) + u2x(k))
            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + u2y(k)           
        enddo
      enddo 
c      print *, 'The stress from wall is done.'

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
cc       trapezoid case
c        do itar = 1,nbodies
c          if (itar .eq. isou) then
c            cycle
c          endif
cc         skip the diagonal term since this was taking care above with
cc         the trapezoid rule with odd even integration
c       
c
cc         loop over target points
c          do k = 1,ninner
cc           loop over source points
c            do j = 1,ninner
c              rx = x((itar-1)*ninner+k) - x((isou-1)*ninner+j)
c              ry = y((itar-1)*ninner+k) - y((isou-1)*ninner+j)
c              rho2 = rx**2.d0 + ry**2.d0
c              rdotn = rx*px((isou-1)*ninner+j) +
c     $                ry*py((isou-1)*ninner+j)
c              routn = rx*py((isou-1)*ninner+j) +
c     $                ry*px((isou-1)*ninner+j)
c              sx = denx(j)
c              sy = deny(j)
c              rdots = rx*sx + ry*sy
c              routs = rx*sy + ry*sx
c
c              E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + 
c     $            5.d-1*
c     $            (2.d0*rdotn*rdots/rho2/rho2 + 
c     $            2.d0*rdots/rho2/rho2*rx*px((isou-1)*ninner+j) +
c     $            2.d0*rdotn/rho2/rho2*rx*sx -
c     $            8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
c     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
c
c              E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) + 
c     $            5.d-1*
c     $            (rdots*routn/rho2/rho2 +
c     $            rdotn*routs/rho2/rho2 -
c     $            8.d0*rdotn*rdots*rx*ry/rho2**3.d0)*
c     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
c
c              E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 
c     $            5.d-1*
c     $            (2.d0*rdotn*rdots/rho2/rho2 + 
c     $            2.d0*rdots/rho2/rho2*ry*py((isou-1)*ninner+j) +
c     $            2.d0*rdotn/rho2/rho2*ry*sy -
c     $            8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
c     $            speed((isou-1)*ninner+j)*twopi/dble(ninner)/pi
c            enddo
c          enddo
c        enddo

c       barycentric case
        do itar = 1,nbodies
          if (itar .eq. isou) then
            cycle
          endif
c         skip the diagonal term since this was taking care above with
c         the trapezoid rule with the correcting liming value at the
c         diagonal

c         loop over target points
          do k = 1,ninner
            xtar(k) = x((itar-1)*ninner+k)
            ytar(k) = y((itar-1)*ninner+k)
            xsou(k) = x((isou-1)*ninner+k)
            ysou(k) = y((isou-1)*ninner+k)            
          enddo
c          print *, itar,isou
            
c           loop over source points
       nbeta = 1
       ninnc = ninner
       nder = 2
       call StokesExteriorDLP(ninner,ninnc,nbeta,xsou,ysou,
     &          denx,deny,px,py,ninner,xtar,ytar,isou,nbodies,nder,
     &          ux,uy,u1x,u1y,u2x,u2y) 
   
c          print *, itar,isou     
          do k = 1,ninner
            E11((itar-1)*ninner+k) = E11((itar-1)*ninner+k) + u1x(k) 
            E12((itar-1)*ninner+k) = E12((itar-1)*ninner+k) +  
     &                               0.5d0*(u1y(k) + u2x(k))
            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + u2y(k)   
          enddo
        enddo

cc       END OF SOURCE POINTS == OBSTACLES, TARGET POINTS ~= OBSTACLE
cc       isou
      enddo
      
c      print *, 'The stress from another bodies is done.'      
c     START OF JUMP ALONG DIAGONAL TERM
      call DCFFTI(ninner,wsave)
      do ibod = 1,nbodies
        do k=1,ninner 
          zden(k) = den(2*nouter+(ibod-1)*2*ninner+k) + 
     $          eye*den(2*nouter+(ibod-1)*2*ninner+k+ninner)
        enddo
        call fourierDiff(ninner,zden)!,wsave)
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
            E22((itar-1)*ninner+k) = E22((itar-1)*ninner+k) + 
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
      subroutine computeShearStress(ninner,nbodies,nouter,x,y,den,ib,
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

      if( ib .eq. 1 ) then
      print *, 'Barycentric deformation'
      call deformation_on_boundary(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)     
      else
      print *, 'Trapezoid deformation'
      call deformation_on_boundary_trap(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
      endif
     
c     compute the deformation tensor on the boundary of the interface

c      if (1 .eq. 1) then
c     
c       open(unit=10,file='output/E11.dat')
c       open(unit=11,file='output/E12.dat')
c       open(unit=12,file='output/E22.dat')       
c       do k = 1,ninner*nbodies
c         write(10,1000) E11(k)
c         write(11,1000) E12(k)
c         write(12,1000) E22(k)         
c       enddo      
c       close(unit=10)
c       close(unit=11)
c       close(unit=12)       
c      endif
c
c

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

 1000 format(E25.16)
 
      end

***********************************************************************
      subroutine computeQoiTargets(ninner,nbodies,nouter,ibary,
     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
c     Compute the velocity, pressure, and vorticity at a set of target 
c     points xtar and ytar.  Target points must be sufficiently far 
c     away from the each boundary since no near-singular integration 
c     is used.  
c     It is just the vanilla trapezoid rule, but will assign a value of
c     zero for points that are inside or close to the boundary
      implicit real*8 (a-h,o-z)

c     input variables
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension xisou(ninner),yisou(ninner)
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
      dimension utar(ntargets),vtar(ntargets),
     &          utmp(ntargets),vtmp(ntargets)
      dimension uxtmp(ntargets),uytmp(ntargets),
     &          vxtmp(ntargets),vytmp(ntargets)     
      dimension press_tar(ntargets)
      dimension vort_tar(ntargets)

c     local variables
      dimension iside(ntargets)
      dimension nnear(nbodies)
      dimension indnear(ntargets,nbodies)
      dimension iup(ntargets,nbodies)
      dimension inear(ntargets)

c     nnear is the number of near points for each body
c     indnear is the index of the points that are close to a 
c     particular body
      call classifyPointsOld(ninner,nbodies,x,y,
     $      ntargets,xtar,ytar,iside,inear)
     
      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
     $    centerx,centery)
c     build inner geometry

      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     build outer geometry

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do j = 1,ntargets
        utar(j) = 0.d0
        vtar(j) = 0.d0
        press_tar(j) = 0.d0
        vort_tar(j) = 0.d0
      enddo

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
          routn = -rx*py0(k) + ry*px0(k)
          routden = -rx*deny(k) + ry*denx(k)
          
          if(ibary .eq. 1) then
             press_tar(j) = press_tar(j) + 
     $         (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $         speed0(k)*twopi/dble(nouter)/pi
          else
            utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $        speed0(k)*twopi/dble(nouter)/pi
            vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $        speed0(k)*twopi/dble(nouter)/pi
            press_tar(j) = press_tar(j) + 
     $        (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $        speed0(k)*twopi/dble(nouter)/pi
            vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*routden)/
     $        rho2**2.d0*speed0(k)*twopi/dble(nouter)/pi
          endif
        enddo
      enddo
      
      if(ibary .eq. 1) then
        nder = 2 
        nbeta = 1
        noutc = nouter
        call StokesInteriorDLP(nouter,noutc,nbeta,xouter,youter,
     &           denx,deny,px0,py0,ntargets,xtar,ytar,nder,
     &           utmp,vtmp,uxtmp,uytmp,vxtmp,vytmp)
           
        do j = 1,ntargets
          utar(j) = utar(j) + utmp(j)       
          vtar(j) = vtar(j) + vtmp(j)
          vort_tar(j) = vort_tar(j) + vxtmp(j)-uytmp(j)
        enddo
      endif 
       
c     Contribution from the density function on each inner obstacle
      do isou = 1,nbodies
        do k = 1,ninner
          denx(k) = den(2*nouter + (isou-1)*2*ninner + k)
          deny(k) = den(2*nouter + (isou-1)*2*ninner + k + ninner)
          xisou(k) = x((isou-1)*ninner + k)
          yisou(k) = y((isou-1)*ninner + k)
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
            routn = -rx*py((isou-1)*ninner+k) + ry*px((isou-1)*ninner+k)
            routden = -rx*deny(k) + ry*denx(k)

            if(ibary .eq. 1) then
              press_tar(j) = press_tar(j) + 
     $          (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $          speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
            else
              utar(j) = utar(j) + rdotn/rho2*rdotden/rho2*rx*
     $          speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
              vtar(j) = vtar(j) + rdotn/rho2*rdotden/rho2*ry*
     $          speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
              press_tar(j) = press_tar(j) + 
     $           (2.d0*rdotn*rdotden/rho2/rho2 - dendotn/rho2)*
     $           speed((isou-1)*ninner+k)*twopi/dble(ninner)/pi
              vort_tar(j) = vort_tar(j) + (rdotden*routn + rdotn*
     $           routden)/rho2**2.d0*speed((isou-1)*ninner+k)*
     $           twopi/dble(ninner)/pi
            endif
          enddo
        enddo
        
        if(ibary .eq. 1) then
         nder = 2 
         nbeta = 1
         ninnc = ninner 
         call StokesExteriorDLP(ninner,ninnc,nbeta,xisou,yisou,
     &         denx,deny,px,py,ntargets,xtar,ytar,isou,nbodies,nder,
     &         utmp,vtmp,uxtmp,uytmp,vxtmp,vytmp)
           
         do j = 1,ntargets
           utar(j) = utar(j) + utmp(j)       
           vtar(j) = vtar(j) + vtmp(j)
           vort_tar(j) = vort_tar(j) + vxtmp(j)-uytmp(j)
         enddo
        endif
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
          vort_tar(j) = vort_tar(j) + (ry*sto1 - rx*sto2)/
     $        rho2/twopi
c         stokeslet contribution

          utar(j) = utar(j) + rot*ry/rho2
          vtar(j) = vtar(j) - rot*rx/rho2
c         rotlet contribution
        enddo
      enddo
c     Add in contribution from Rotlets and Stokeslets

      do k = 1,ntargets
c        if (iside(k) .eq. 0 .or. inear(k) .eq. 1) then
        if (iside(k) .eq. 0) then
          utar(k) = 0.d0
          vtar(k) = 0.d0
          press_tar(k) = 0.d0
          vort_tar(k) = 0.d0
        endif
      enddo

c      do j=1,nbodies
c        do k = 1,nnear(j)
c          if (iside(indnear(k,j)) .ne. 0) then
c            utar(indnear(k,j)) = dble(iup(k,j))
c            vtar(indnear(k,j)) = dble(iup(k,j))
c          endif
c        enddo
c      enddo



      end
      
c***********************************************************************
c      subroutine computeQoiTargetsTrap(ninner,nbodies,nouter,
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
c      dimension inear(ntargets)
c
cc     nnear is the number of near points for each body
cc     indnear is the index of the points that are close to a 
cc     particular body
c      call classifyPointsOld(ninner,nbodies,x,y,
c     $      ntargets,xtar,ytar,iside,inear)
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
cc      do k = 1,ntargets
ccc        if (iside(k) .eq. 0 .or. inear(k) .eq. 1) then
cc        if (iside(k) .eq. 0) then
cc          utar(k) = 0.d0
cc          vtar(k) = 0.d0
cc          press_tar(k) = 0.d0
cc          vort_tar(k) = 0.d0
cc        endif
cc      enddo
c
c
c
cc      do j=1,nbodies
cc        do k = 1,nnear(j)
cc          if (iside(indnear(k,j)) .ne. 0) then
cc            utar(indnear(k,j)) = dble(iup(k,j))
cc            vtar(indnear(k,j)) = dble(iup(k,j))
cc          endif
cc        enddo
cc      enddo
cc
c
c
c      end

c***********************************************************************
c      subroutine computeQoiTargets(ninner,nbodies,nouter,
c     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
cc     Compute the velocity, pressure, and vorticity at a set of target 
cc     points xtar and ytar.
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
c      dimension iup(ntargets)
c      dimension xtar_level(ntargets),ytar_level(ntargets)
c      dimension utar_level(ntargets),vtar_level(ntargets)
c      dimension press_tar_level(ntargets),vort_tar_level(ntargets)
c      dimension xup(16*max(ninner,nouter))
c      dimension yup(16*max(ninner,nouter))
c      dimension pxup(16*max(ninner,nouter))
c      dimension pyup(16*max(ninner,nouter))
c      dimension curup(16*max(ninner,nouter))
c      dimension speedup(16*max(ninner,nouter))
c      dimension denxup(16*max(ninner,nouter))
c      dimension denyup(16*max(ninner,nouter))
c      complex *16 eye
c      complex *16 zin(16*max(ninner,nouter))
c      complex *16 zout(16*max(ninner,nouter))
c
c      pi = 4.d0*datan(1.d0)
c      twopi = 2.d0*pi
c      eye = (0.d0,1.d0)
c
c      do j = 1,ntargets
c        utar(j) = 0.d0
c        vtar(j) = 0.d0
c        press_tar(j) = 0.d0
c        vort_tar(j) = 0.d0
c      enddo
c
c      call inner_geometry(ninner,nbodies,x,y,x,y,px,py,cur,speed,
c     $    centerx,centery)
cc     build inner geometry
c
cc     build inner geometry
c      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
cc     build outer geometry
c
cc     Contribution from the density function on the solid wall
c      do k=1,nouter
c        denx(k) = den(k)
c        deny(k) = den(k+nouter)
c      enddo
cc     Density function defined on the outer geometry
c
c
c      do j = 1,nouter
c        xup(j) = xouter(j)
c        yup(j) = youter(j)
c        denxup(j) = denx(j)
c        denyup(j) = deny(j)
c      enddo
c
c      call classifyPoints(nouter,xup,yup,
c     $    ntargets,xtar,ytar,iside,nnear,iup)
c
c      do k = 1,ntargets
c        if (iside(k) .eq. 1 .or. iup(k) .gt. 16) then
c          utar(k) = 1.0d20
c          vtar(k) = 1.0d20
c          press_tar(k) = 1.0d20
c          vort_tar(k) = 1.0d20
c        endif
c      enddo
c
cc     START OF DOING POINTS THAT REQUIRE NO UPSAMPLING
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 1 .and. iside(j) .eq. 0) then
c          xtar_level(icount) = xtar(j)
c          ytar_level(icount) = ytar(j)
c          icount = icount + 1
c        endif
c      enddo
c
c      call evalLPouter(nouter,xup,yup,
c     $  denxup,denyup,icount-1,xtar_level,ytar_level,
c     $  utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 1 .and. iside(j) .eq. 0) then
c          utar(j) = utar(j) + utar_level(icount)
c          vtar(j) = vtar(j) + vtar_level(icount)
c          press_tar(j) = press_tar(j) + press_tar_level(icount)
c          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c          icount = icount + 1
c        endif
c      enddo
cc     END OF DOING POINTS THAT REQUIRE NO UPSAMPLING
c
cc     START OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 2 .and. iside(j) .eq. 0) then
c          xtar_level(icount) = xtar(j)
c          ytar_level(icount) = ytar(j)
c          icount = icount + 1
c        endif
c      enddo
c
c      do k = 1,nouter
c        zin(k) = xup(k) + eye*yup(k)
c      enddo
c      call fourierUpsample(nouter,2,zin,zout)
c      do k = 1,2*nouter
c        xup(k) = dreal(zout(k))
c        yup(k) = dimag(zout(k))
c      enddo
cc     upsample geometry by a factor of 2
c      do k = 1,nouter
c        zin(k) = denxup(k) + eye*denyup(k)
c      enddo
c      call fourierUpsample(nouter,2,zin,zout)
c      do k = 1,2*nouter
c        denxup(k) = dreal(zout(k))
c        denyup(k) = dimag(zout(k))
c      enddo
cc     upsample density by a factor of 2
c
c      call evalLPouter(2*nouter,xup,yup,
c     $  denxup,denyup,icount-1,xtar_level,ytar_level,
c     $  utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 2 .and. iside(j) .eq. 0) then
c          utar(j) = utar(j) + utar_level(icount)
c          vtar(j) = vtar(j) + vtar_level(icount)
c          press_tar(j) = press_tar(j) + press_tar_level(icount)
c          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c          icount = icount + 1
c        endif
c      enddo
cc     END OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
c
c
cc     START OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 4 .and. iside(j) .eq. 0) then
c          xtar_level(icount) = xtar(j)
c          ytar_level(icount) = ytar(j)
c          icount = icount + 1
c        endif
c      enddo
c
c      do k = 1,2*nouter
c        zin(k) = xup(k) + eye*yup(k)
c      enddo
c      call fourierUpsample(2*nouter,2,zin,zout)
c      do k = 1,4*nouter
c        xup(k) = dreal(zout(k))
c        yup(k) = dimag(zout(k))
c      enddo
cc     upsample geometry by a factor of 2
c      do k = 1,2*nouter
c        zin(k) = denxup(k) + eye*denyup(k)
c      enddo
c      call fourierUpsample(2*nouter,2,zin,zout)
c      do k = 1,4*nouter
c        denxup(k) = dreal(zout(k))
c        denyup(k) = dimag(zout(k))
c      enddo
cc     upsample density by a factor of 2
c
c      call evalLPouter(4*nouter,xup,yup,
c     $  denxup,denyup,icount-1,xtar_level,ytar_level,
c     $  utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 4 .and. iside(j) .eq. 0) then
c          utar(j) = utar(j) + utar_level(icount)
c          vtar(j) = vtar(j) + vtar_level(icount)
c          press_tar(j) = press_tar(j) + press_tar_level(icount)
c          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c          icount = icount + 1
c        endif
c      enddo
cc     END OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
c
c
cc     START OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 8 .and. iside(j) .eq. 0) then
c          xtar_level(icount) = xtar(j)
c          ytar_level(icount) = ytar(j)
c          icount = icount + 1
c        endif
c      enddo
c
c      do k = 1,4*nouter
c        zin(k) = xup(k) + eye*yup(k)
c      enddo
c      call fourierUpsample(4*nouter,2,zin,zout)
c      do k = 1,8*nouter
c        xup(k) = dreal(zout(k))
c        yup(k) = dimag(zout(k))
c      enddo
cc     upsample geometry by a factor of 2
c      do k = 1,4*nouter
c        zin(k) = denxup(k) + eye*denyup(k)
c      enddo
c      call fourierUpsample(4*nouter,2,zin,zout)
c      do k = 1,8*nouter
c        denxup(k) = dreal(zout(k))
c        denyup(k) = dimag(zout(k))
c      enddo
cc     upsample density by a factor of 2
c
c      call evalLPouter(8*nouter,xup,yup,
c     $  denxup,denyup,icount-1,xtar_level,ytar_level,
c     $  utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 8 .and. iside(j) .eq. 0) then
c          utar(j) = utar(j) + utar_level(icount)
c          vtar(j) = vtar(j) + vtar_level(icount)
c          press_tar(j) = press_tar(j) + press_tar_level(icount)
c          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c          icount = icount + 1
c        endif
c      enddo
cc     END OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
c
cc     START OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 16 .and. iside(j) .eq. 0) then
c          xtar_level(icount) = xtar(j)
c          ytar_level(icount) = ytar(j)
c          icount = icount + 1
c        endif
c      enddo
c
c      do k = 1,8*nouter
c        zin(k) = xup(k) + eye*yup(k)
c      enddo
c      call fourierUpsample(8*nouter,2,zin,zout)
c      do k = 1,16*nouter
c        xup(k) = dreal(zout(k))
c        yup(k) = dimag(zout(k))
c      enddo
cc     upsample geometry by a factor of 2
c      do k = 1,8*nouter
c        zin(k) = denxup(k) + eye*denyup(k)
c      enddo
c      call fourierUpsample(8*nouter,2,zin,zout)
c      do k = 1,16*nouter
c        denxup(k) = dreal(zout(k))
c        denyup(k) = dimag(zout(k))
c      enddo
cc     upsample density by a factor of 2
c
c      call evalLPouter(16*nouter,xup,yup,
c     $  denxup,denyup,icount-1,xtar_level,ytar_level,
c     $  utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c      icount = 1
c      do j=1,ntargets
c        if (iup(j) .eq. 16 .and. iside(j) .eq. 0) then
c          utar(j) = utar(j) + utar_level(icount)
c          vtar(j) = vtar(j) + vtar_level(icount)
c          press_tar(j) = press_tar(j) + press_tar_level(icount)
c          vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c          icount = icount + 1
c        endif
c      enddo
cc     END OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
c
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
c      do ibod = 1,nbodies
c        do k = 1,ninner
c          denx(k) = den(2*nouter + (ibod-1)*2*ninner + k)
c          deny(k) = den(2*nouter + (ibod-1)*2*ninner + k + ninner)
c        enddo
c
c        do j = 1,ninner
c          xup(j) = x((ibod-1)*ninner + j)
c          yup(j) = y((ibod-1)*ninner + j)
c          denxup(j) = denx(j)
c          denyup(j) = deny(j)
c        enddo
cc
c        call classifyPoints(ninner,xup,yup,
c     $    ntargets,xtar,ytar,iside,nnear,iup)
c
cc       START OF DOING POINTS THAT REQUIRE NO UPSAMPLING
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 1) then
c            xtar_level(icount) = xtar(j)
c            ytar_level(icount) = ytar(j)
c            icount = icount + 1
c          endif
c        enddo
c
c        call evalLPbodies(ninner,xup,yup,
c     $    denxup,denyup,icount-1,xtar_level,ytar_level,
c     $    utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 1) then
c            utar(j) = utar(j) + utar_level(icount)
c            vtar(j) = vtar(j) + vtar_level(icount)
c            press_tar(j) = press_tar(j) + press_tar_level(icount)
c            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c            icount = icount + 1
c          endif
c        enddo
cc       END OF DOING POINTS THAT REQUIRE NO UPSAMPLING
c
cc       START OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 2) then
c            xtar_level(icount) = xtar(j)
c            ytar_level(icount) = ytar(j)
c            icount = icount + 1
c          endif
c        enddo
c
c        do k = 1,ninner
c          zin(k) = xup(k) + eye*yup(k)
c        enddo
c        call fourierUpsample(ninner,2,zin,zout)
c        do k = 1,2*ninner
c          xup(k) = dreal(zout(k))
c          yup(k) = dimag(zout(k))
c        enddo
cc       upsample geometry by a factor of 2
c        do k = 1,ninner
c          zin(k) = denxup(k) + eye*denyup(k)
c        enddo
c        call fourierUpsample(ninner,2,zin,zout)
c        do k = 1,2*ninner
c          denxup(k) = dreal(zout(k))
c          denyup(k) = dimag(zout(k))
c        enddo
cc       upsample density by a factor of 2
c
c        call evalLPbodies(2*ninner,xup,yup,
c     $    denxup,denyup,icount-1,xtar_level,ytar_level,
c     $    utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 2) then
c            utar(j) = utar(j) + utar_level(icount)
c            vtar(j) = vtar(j) + vtar_level(icount)
c            press_tar(j) = press_tar(j) + press_tar_level(icount)
c            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c            icount = icount + 1
c          endif
c        enddo
cc       END OF DOING POINTS THAT REQUIRE 2X UPSAMPLING
c
cc       START OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 4) then
c            xtar_level(icount) = xtar(j)
c            ytar_level(icount) = ytar(j)
c            icount = icount + 1
c          endif
c        enddo
c
c        do k = 1,2*ninner
c          zin(k) = xup(k) + eye*yup(k)
c        enddo
c        call fourierUpsample(2*ninner,2,zin,zout)
c        do k = 1,4*ninner
c          xup(k) = dreal(zout(k))
c          yup(k) = dimag(zout(k))
c        enddo
cc       upsample geometry by a factor of 2
c        do k = 1,2*ninner
c          zin(k) = denxup(k) + eye*denyup(k)
c        enddo
c        call fourierUpsample(2*ninner,2,zin,zout)
c        do k = 1,4*ninner
c          denxup(k) = dreal(zout(k))
c          denyup(k) = dimag(zout(k))
c        enddo
cc       upsample density by a factor of 2
c
c        call evalLPbodies(4*ninner,xup,yup,
c     $    denxup,denyup,icount-1,xtar_level,ytar_level,
c     $    utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 4) then
c            utar(j) = utar(j) + utar_level(icount)
c            vtar(j) = vtar(j) + vtar_level(icount)
c            press_tar(j) = press_tar(j) + press_tar_level(icount)
c            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c            icount = icount + 1
c          endif
c        enddo
cc       END OF DOING POINTS THAT REQUIRE 4X UPSAMPLING
c
cc       START OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 8) then
c            xtar_level(icount) = xtar(j)
c            ytar_level(icount) = ytar(j)
c            icount = icount + 1
c          endif
c        enddo
c
c        do k = 1,4*ninner
c          zin(k) = xup(k) + eye*yup(k)
c        enddo
c        call fourierUpsample(4*ninner,2,zin,zout)
c        do k = 1,8*ninner
c          xup(k) = dreal(zout(k))
c          yup(k) = dimag(zout(k))
c        enddo
cc       upsample geometry by a factor of 2
c        do k = 1,4*ninner
c          zin(k) = denxup(k) + eye*denyup(k)
c        enddo
c        call fourierUpsample(4*ninner,2,zin,zout)
c        do k = 1,8*ninner
c          denxup(k) = dreal(zout(k))
c          denyup(k) = dimag(zout(k))
c        enddo
cc       upsample density by a factor of 2
c
c        call evalLPbodies(8*ninner,xup,yup,
c     $    denxup,denyup,icount-1,xtar_level,ytar_level,
c     $    utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 8) then
c            utar(j) = utar(j) + utar_level(icount)
c            vtar(j) = vtar(j) + vtar_level(icount)
c            press_tar(j) = press_tar(j) + press_tar_level(icount)
c            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c            icount = icount + 1
c          endif
c        enddo
cc       END OF DOING POINTS THAT REQUIRE 8X UPSAMPLING
c
cc       START OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 16) then
c            xtar_level(icount) = xtar(j)
c            ytar_level(icount) = ytar(j)
c            icount = icount + 1
c          endif
c        enddo
c
c        do k = 1,8*ninner
c          zin(k) = xup(k) + eye*yup(k)
c        enddo
c        call fourierUpsample(8*ninner,2,zin,zout)
c        do k = 1,16*ninner
c          xup(k) = dreal(zout(k))
c          yup(k) = dimag(zout(k))
c        enddo
cc       upsample geometry by a factor of 2
c        do k = 1,8*ninner
c          zin(k) = denxup(k) + eye*denyup(k)
c        enddo
c        call fourierUpsample(8*ninner,2,zin,zout)
c        do k = 1,16*ninner
c          denxup(k) = dreal(zout(k))
c          denyup(k) = dimag(zout(k))
c        enddo
cc       upsample density by a factor of 2
c
c        call evalLPbodies(16*ninner,xup,yup,
c     $    denxup,denyup,icount-1,xtar_level,ytar_level,
c     $    utar_level,vtar_level,press_tar_level,vort_tar_level)
c
c        icount = 1
c        do j=1,ntargets
c          if (iup(j) .eq. 16) then
c            utar(j) = utar(j) + utar_level(icount)
c            vtar(j) = vtar(j) + vtar_level(icount)
c            press_tar(j) = press_tar(j) + press_tar_level(icount)
c            vort_tar(j) = vort_tar(j) + vort_tar_level(icount)
c            icount = icount + 1
c          endif
c        enddo
cc       END OF DOING POINTS THAT REQUIRE 16X UPSAMPLING
c
c        do j=1,ntargets
c          if (iside(j) .eq. 0 .or. iup(j) .gt. 16) then
c            utar(j) = 1.0d20
c            vtar(j) = 1.0d20
c            press_tar(j) = 1.0d20
c            vort_tar(j) = 1.0d20
c          endif
c        enddo
c      enddo
c
c      do k = 1,ntargets
c        if (abs(utar(k)) .ge. 1.0d8) then
c          utar(k) = 0.d0
c          vtar(k) = 0.d0
c          press_tar(k) = 0.d0
c          vort_tar(k) = 0.d0
c        endif
c      enddo
c
c
c      end
c
c***********************************************************************
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
        call fourierDiff(ninner,zden)!,wsave)
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

c*******************************************************************************8      
      subroutine StokesInteriorDLP(npts,nptsc,nbeta,xsou,ysou,
     &          dens1,dens2,nx,ny,ntar,xtar,ytar,nder,
     &          utar1x,utar1y,u1xtar,u1ytar,u2xtar,u2ytar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts),
     &          xsouc(nptsc),ysouc(nptsc)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     &            ddutar(ntar)
c      real *8 uxtar(ntar),uytar(ntar)
      real *8 utar1x(ntar),utar1y(ntar)      
      real *8 u1xtar(ntar),u1ytar(ntar),
     &        u2xtar(ntar),u2ytar(ntar)  
      real*8 nx(npts),ny(npts)
      complex *16 zsou(npts), outnor



      dimension dens1(npts), dens2(npts), sdotd(nptsc),
     &          dens1c(nptsc), dens2c(nptsc)
      complex *16 tau1(npts), tau2(npts)
      complex *16 eye, czero



      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.d0)
      czero = (0.d0, 0.d0)

c       DLP density function
c        Input data

c       The target points
c        Input data


c      a coarse group of collocation points and the density at them      

       Do j=1, nptsc
          k=nbeta*(j-1)+1
          xsouc(j) = xsou(k)
          ysouc(j) = ysou(k)
          dens1c(j) = dens1(k)
          dens2c(j) = dens2(k)
          sdotd(j) = xsouc(j)*dens1c(j) + ysouc(j)*dens2c(j)
       ENDDO

 
c      Find the outward normal vectors on the boundary
c        Input data nx,ny

c      Find the functions for the first part of integrations in equation (2.9)
       do j=1,npts
         outnor = czero
         outnor = nx(j)+ eye*ny(j)
         tau1(j) = nx(j)*(dens1(j) + eye*dens2(j))/outnor
         tau2(j) = ny(j)*(dens1(j) + eye*dens2(j))/outnor
       enddo


c      Find DLP by doing four integrations in equation (2.9) of
c      Barnett, Wu, Veerapaneni

c      The Second, Third and, Forth parts
       

      do j = 1,ntar
        utar1x(j) = 0.d0
        utar1y(j) = 0.d0
        u1xtar(j) = 0.d0
        u1ytar(j) = 0.d0
        u2xtar(j) = 0.d0
        u2ytar(j) = 0.d0       
      enddo

      call compute_bdval_in(nptsc,xsou,ysou,dcmplx(sdotd),bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)
     
      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(dutar(j))
        utar1y(j) = utar1y(j) - dimag(dutar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) - dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) - dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) - dreal(ddutar(j))
        endif
      enddo     
     
     
      call compute_bdval_in(nptsc,xsou,ysou,dcmplx(dens1c),bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)     

      do j = 1,ntar
        utar1x(j) = utar1x(j) - xtar(j)*dreal(dutar(j))
        utar1y(j) = utar1y(j) + xtar(j)*dimag(dutar(j))
        
        if(nder .eq. 2) then
        u1xtar(j) = u1xtar(j) - dreal(dutar(j)) - xtar(j)*
     $              dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) + xtar(j)*dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) + dimag(dutar(j)) + xtar(j)*
     $              dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) + xtar(j)*dreal(ddutar(j))
        endif
        
      enddo        

      call compute_bdval_in(nptsc,xsou,ysou,dcmplx(dens2c),bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)     

      do j = 1,ntar
        utar1x(j) = utar1x(j) - ytar(j)*dreal(dutar(j))
        utar1y(j) = utar1y(j) + ytar(j)*dimag(dutar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) - ytar(j)*dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) - dreal(dutar(j)) + ytar(j)*
     $              dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) + ytar(j)*dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) + dimag(dutar(j)) + ytar(j)*
     $              dreal(ddutar(j))
        endif
      enddo         

      call compute_bdval_in(npts,xsou,ysou,tau1,bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)

      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(dutar(j))
        u1ytar(j) = u1ytar(j) - dimag(dutar(j))
        endif
      enddo

      call compute_bdval_in(npts,xsou,ysou,tau2,bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)

      do j = 1,ntar
        utar1y(j) = utar1y(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u2xtar(j) = u2xtar(j) + dreal(dutar(j))
        u2ytar(j) = u2ytar(j) - dimag(dutar(j))
        endif
      enddo
      
      end
      
******************************************************************
      subroutine compute_bdval_in(npts,xsou,ysou,density,bdval)
      implicit real*8 (a-h,o-z)
      
      dimension xsou(npts),ysou(npts)
      complex *16 density(npts)
      complex *16 bdval(npts)

      complex *16 zsou(npts),dzsou(npts)
      complex *16 eye,czero
      complex *16 ddensity(npts)

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.d0)
      czero = (0.d0,0.d0)

      do j = 1,npts
        ddensity(j) = density(j)
      enddo
      call fourierDiff(npts,ddensity)
c     compute derivative of the density function

      do k=1,npts
        zsou(k) = xsou(k) + eye*ysou(k)
        dzsou(k) = zsou(k)*dtheta
      enddo
c     put geometry in complex variables
      call fourierDiff(npts,dzsou)

      do j=1,npts
        bdval(j) = czero
        do k=1,npts
          if (k .ne. j) then
            bdval(j) = bdval(j)-dcmplx(density(k) - density(j))/
     $          (zsou(k) - zsou(j))*dzsou(k)
          else
            bdval(j) = bdval(j) - twopi/dble(npts)*ddensity(j)
          endif
        enddo
        bdval(j) = -eye*bdval(j)/twopi - density(j)
c        bdval(j) = -bdval(j)
      enddo

      return
      end

******************************************************************
      subroutine laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $      ntar,xtar,ytar,utar,nder,dutar,ddutar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     $             ddutar(ntar)

      complex *16 eye,czero
      complex *16 zsou(npts),dzsou(npts)
      complex *16 den1(npts),den2(npts)
      complex *16 num,den,num1,den3
      complex *16 ztar
      complex *16 ftar1,dftar1
      complex *16 ftar2,dftar2

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.0d0)
      czero = (0.d0,0.0d0)

      dtheta = twopi/dble(npts)
      do k=1,npts
        zsou(k) = xsou(k) + eye*ysou(k)
        dzsou(k) = zsou(k)*dtheta
      enddo
c     put geometry in complex variables
      call fourierDiff(npts,dzsou)

c     FMM for the denominator
c     FMM for the first numerator ... also compute the gradient
c     FMM for the second numerator with weight function bdval .. also
compute the gradient which you need for the third part
c     for j=1,ntar
c       num(j) = num(j) + utar(j)*grad(j)
c     enddo

c     den2 = FMMcall()

      do j = 1,ntar
        ztar = xtar(j) + eye*ytar(j)
        num = czero
        den = czero
        do k = 1,npts
          den1(k) = zsou(k) - ztar
          num = num + bdval(k)/den1(k)*dzsou(k)
          den = den + 1.d0/den1(k)*dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.5) of
c       of Barnett, Wu, Veerapaneni

        utar(j) = num/den

c       Start computing derivative (equation (3.10)).  Term in
c       denominator is unchanged
        num = czero
        do k = 1,npts
          den2(k) = den1(k)*den1(k)          
          num = num + (bdval(k) - utar(j))/
     $          den2(k)*dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.5) of
c       Barnett, Wu, Veerapaneni
        dutar(j) = num/den

c       Computing the second derivative of DLP for the stress tensor
        if(nder .eq. 2) then
          num = czero
          num1 = czero
          den3 = czero
          do k=1,npts
            den3 = den2(k)*den1(k)
            num1 = num1+(bdval(k)-utar(j)-dutar(j)*(zsou(k)-ztar))/
     $             den3*dzsou(k)
          enddo
          ddutar(j)=2.d0*num1/den
        endif
      enddo

c      do k = 1,ntar
c        print*,den(k) - den2(k)
c      enddo

      end

c****************************************************************
      subroutine StokesExteriorDLP(npts,nptsc,nbeta,xsou,ysou,
     &          densx,densy,px,py,ntar,xtar,ytar,isou,n,nder,
     &          utar1x,utar1y,u1xtar,u1ytar,u2xtar,u2ytar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts),
     &          xsouc(nptsc),ysouc(nptsc)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     &            ddutar(ntar)

      dimension utar1x(ntar),utar1y(ntar)
      dimension u1xtar(ntar),u1ytar(ntar),
     &          u2xtar(ntar),u2ytar(ntar)      

      dimension px(npts*n),py(npts*n),px1(npts),py1(npts)

      complex *16 zsou(npts), outnor

c      real*8 omega
      dimension densx(npts), densy(npts), sdotd(npts)!,
!     &          dens1c(nptsc), dens2c(nptsc)
      complex *16 tau1(npts), tau2(npts)
      complex *16 eye, czero



c     the center of the exterior boundary when doing parametrization of boundary 
c      xc=0.d0
c      yc=0.d0

      twopi = 8.d0*datan(1.d0)
      pi = 4.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.d0)
      czero = (0.d0, 0.d0)


      do j=1,npts
        outnor = czero
        px1(j) = px((isou-1)*npts+j)
        py1(j) = py((isou-1)*npts+j)
        outnor = px1(j) + eye*py1(j)
        tau1(j) = (densx(j) + eye*densy(j))*px1(j)/outnor
        tau2(j) = (densx(j) + eye*densy(j))*py1(j)/outnor
        sdotd(j) = xsou(j)*densx(j) + ysou(j)*densy(j)         
      enddo      


c      Find DLP by doing four integrations in equation (2.9) of
c      Barnett, Wu, Veerapaneni

c      The Second, Third and, Forth parts
      do j = 1,ntar
        utar1x(j) = 0.d0
        utar1y(j) = 0.d0
        u1xtar(j) = 0.d0
        u1ytar(j) = 0.d0
        u2xtar(j) = 0.d0
        u2ytar(j) = 0.d0         
      enddo

      call compute_bdval_ex(npts,xsou,ysou,dcmplx(sdotd),bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)
     
      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(dutar(j))
        utar1y(j) = utar1y(j) - dimag(dutar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) - dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) - dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) - dreal(ddutar(j))
        endif
      enddo     
     
     
      call compute_bdval_ex(npts,xsou,ysou,dcmplx(densx),bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)     

      do j = 1,ntar
        utar1x(j) = utar1x(j) - xtar(j)*dreal(dutar(j))
        utar1y(j) = utar1y(j) + xtar(j)*dimag(dutar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) - dreal(dutar(j)) - xtar(j)*
     $              dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) + xtar(j)*dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) + dimag(dutar(j)) + xtar(j)*
     $              dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) + xtar(j)*dreal(ddutar(j))
        endif
        
      enddo        

      call compute_bdval_ex(npts,xsou,ysou,dcmplx(densy),bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)     

      do j = 1,ntar
        utar1x(j) = utar1x(j) - ytar(j)*dreal(dutar(j))
        utar1y(j) = utar1y(j) + ytar(j)*dimag(dutar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) - ytar(j)*dreal(ddutar(j))
        u1ytar(j) = u1ytar(j) - dreal(dutar(j)) + ytar(j)*
     $              dimag(ddutar(j))
        u2xtar(j) = u2xtar(j) + ytar(j)*dimag(ddutar(j))
        u2ytar(j) = u2ytar(j) + dimag(dutar(j)) + ytar(j)*
     $              dreal(ddutar(j))
        endif
        
      enddo         
      call compute_bdval_ex(npts,xsou,ysou,tau1,bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)
      
      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(dutar(j))
        u1ytar(j) = u1ytar(j) - dimag(dutar(j))
        endif
        
      enddo

      call compute_bdval_ex(npts,xsou,ysou,tau2,bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,nder,dutar,ddutar)
      
      do j = 1,ntar
        utar1y(j) = utar1y(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u2xtar(j) = u2xtar(j) + dreal(dutar(j))
        u2ytar(j) = u2ytar(j) - dimag(dutar(j))
        endif
        
      enddo 

      return
      end

c******************************************************************
      subroutine compute_bdval_ex(npts,xsou,ysou,density,bdval)
      implicit real*8 (a-h,o-z)
      
      dimension xsou(npts),ysou(npts)
      complex *16 density(npts)
      complex *16 bdval(npts)

      complex *16 zsou(npts),dzsou(npts)
      complex *16 dummy_var,eye,czero
      complex *16 ddensity(npts)

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.d0)
      czero = (0.d0,0.d0)

      do j = 1,npts
        ddensity(j) = dcmplx(density(j))
      enddo
      call fourierDiff(npts,ddensity)
c     compute derivative of the density function

      do k=1,npts
        zsou(k) = xsou(k) + eye*ysou(k)
        dzsou(k) = zsou(k)*dtheta
      enddo
c     put geometry in complex variables
      call fourierDiff(npts,dzsou)

      do j=1,npts
        bdval(j) = czero
        do k=1,npts
          if (k .ne. j) then
            bdval(j) = bdval(j)-dcmplx(density(k) - density(j))/
     $          (zsou(k) - zsou(j))*dzsou(k)
          else
            bdval(j) = bdval(j) - twopi/dble(npts)*ddensity(j)
          endif
        enddo
        bdval(j) = -eye*bdval(j)/twopi 
        bdval(j) = -bdval(j) 
c     The formula (4.3) for exterior limit v^+
      enddo

      return
      end





c******************************************************************

      subroutine laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $      ntar,xtar,ytar,utar,nder,dutar,ddutar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     &            ddutar(ntar) 

      complex *16 eye,czero,a
      complex *16 zsou(npts),dzsou(npts)
      complex *16 den1(npts), den2(npts)
      complex *16 num,den,num1,den3
      complex *16 ztar
      complex *16 ftar1,dftar1
      complex *16 ftar2,dftar2

      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(npts)
      eye = (0.d0,1.d0)
      czero = (0.d0,0.d0)
        x=0.d0
        y=0.d0
      do k=1,npts
        x = x + xsou(k)/dble(npts)
        y = y + ysou(k)/dble(npts)
      enddo
      
      a = x + eye*y
c      print *, x,y,a
      dtheta = twopi/dble(npts)
      do k=1,npts
        zsou(k) = xsou(k) + eye*ysou(k)
        dzsou(k) = zsou(k)*dtheta
      enddo
c     put geometry in complex variables
      call fourierDiff(npts,dzsou)

      do j = 1,ntar
        ztar = xtar(j) + eye*ytar(j)
        num = czero
        den = czero
        do k = 1,npts
          den1(k) = czero
          den1(k) = (zsou(k) - ztar)
          num = num + bdval(k)/den1(k)*dzsou(k)
c          den = den + 1.d0/((zsou(k) - dcmplx(a))
c     $          *(zsou(k) - ztar))*dzsou(k)
          den = den + 1.d0/den1(k)*dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.8) of
c       of Barnett, Wu, Veerapaneni and let a=(0,0) 

c        utar(j) = num/(ztar - dcmplx(a))*den
        utar(j) = num/(-twopi*eye+den)

c       Start computing derivative (equation (3.12)).  Term in
c       denominator is unchanged
        num = czero
        den = czero
        tmp1 = czero
        tmp2 = czero
        do k = 1,npts
          den2(k) = czero
          den2(k) = den1(k)*den1(k)         
        if( CDABS(zsou(k)-ztar) .gt. 1.d-10) then
          num = num + (bdval(k) - utar(j))/
     $          den2(k)*dzsou(k)
          else
c          print *,1
          do l=1,npts
             if( l .ne. k) then
             tmp1 = tmp1 + (bdval(k)*(zsou(k)-dcmplx(a))/
     $          (zsou(l)-dcmplx(a))-bdval(l))*dzsou(l)/(zsou(l)-ztar)
             endif
             tmp2 = tmp2 + 1.d0/((zsou(l) - dcmplx(a))
     $          *(zsou(l) - ztar))*dzsou(l)           
           enddo   
           num = num + (tmp1/tmp2 - (zsou(k) - ztar)*bdval(k))/
     $          (zsou(k) - ztar)**2.d0/(ztar-dcmplx(a))*dzsou(k)
           endif
           den = den +  1.d0/((zsou(k) - dcmplx(a))*(zsou(k) - ztar))
     $          *dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.8) of
c       Barnett, Wu, Veerapaneni
        dutar(j) = num/(ztar - dcmplx(a))/den
        
c       Computing the second derivative of DLP for the stress tensor
        if(nder .eq. 2) then
        num = czero
        num1 = czero
          Do k=1,npts
            den3=den2(k)*den1(k)
            num1 = num1+(bdval(k)-utar(j)-dutar(j)*(zsou(k)-ztar))/
     $               den3*dzsou(k)
          enddo
          ddutar(j)=2.d0*num1/(ztar - dcmplx(a))/den
        endif
      enddo

      return
      end
     
