      subroutine stokesSolver(nninner,nnbodies,nntargets,ifmm,xx,yy,
     $    xxtar,yytar,utar,vtar,press_tar,shear_stress)
c     Input x and y coordinates and return the shear_stress on the
c     boundary.  Outer wall is used to enclose the inner boundary so
c     that Stokes paradox is avoided
      implicit real*8 (a-h,o-z)

      dimension xx(nninner*nnbodies),yy(nninner*nnbodies)
c     x and y coordinates of obstacle
      dimension xxtar(nntargets),yytar(nntargets)
c     x and y coordinates of target locations where velocity and
c     pressure need to be evaluted

      parameter (nmax = 2**15)
      parameter (maxbodies = 50)
c     max points on the boundary of the obstacle      
      parameter (ntarmax = 20000)
      parameter (maxl = 3000, liwork = 30)
      parameter (lrwork = 10 + nmax*(maxl+6) + 
     $     maxl*(maxl+3))
c     maximum size of workspaces for GMRES

      dimension x(nmax),y(nmax)
      dimension centerx(maxbodies),centery(maxbodies)
      dimension px(nmax),py(nmax)
c     x and y coordinates of the normal of the obstacle
      dimension cur(nmax),speed(nmax)
c     Jacobian and curvature of the geometry

      dimension xtar(ntarmax),ytar(ntarmax)
c     Target locations
      dimension utar(nntargets),vtar(nntargets)
      dimension press_tar(nntargets)
c     Velocity and pressure at target locations

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
c     x and y coordinates of the density function
      dimension E11(nninner*nnbodies)
      dimension E12(nninner*nnbodies)
      dimension E22(nninner*nnbodies)
c     components of the deformation gradient on the obstacle
      dimension shear_stress(nninner*nnbodies)
c     shear stress on the obstacle

      common /geometry/x,y,centerx,centery,px,py,cur,speed,
     $    ninner,nbodies
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter
c     global variables for the geometry which are needed by the external
c     matvec routine in gmres

      ninner = nninner
      nbodies = nnbodies
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
     $      gmwork,lrwork,igwork,liwork,maxl,ifmm)
c     solve for the density function with GMRES

c      open(unit=1,file='output/den.dat')
c      do k = 1,2*nouter+2*ninner*nbodies+3*nbodies
c        write(1,1000) den(k)
c      enddo
c      close(1)

      call deformation_on_boundary(ninner,nbodies,x,y,
     $    centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor on the boundary of the interface

      call compute_shear_stress(ninner,nbodies,px,py,
     $    E11,E12,E22,shear_stress)
c     Use the deformation tensor to compute the shear stress

      ntargets = nntargets
      call eval_velocity(ninner,nbodies,x,y,centerx,centery,
     $ px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $ ntargets,xtar,ytar,xxtar,yytar,utar,vtar,press_tar)
c     evaluate velocity and pressure at target points

c      open(unit=1,file='output/xtar.dat')
c      open(unit=2,file='output/ytar.dat')
c      open(unit=3,file='output/utar.dat')
c      open(unit=4,file='output/vtar.dat')
c      open(unit=5,file='output/press_tar.dat')
c      do k = 1,ntargets
c        write(1,1000) xtar(k)
c        write(2,1000) ytar(k)
c        write(3,1000) utar(k)
c        write(4,1000) vtar(k)
c        write(5,1000) press_tar(k)
c      enddo
c      close(unit=1)
c      close(unit=2)
c      close(unit=3)
c      close(unit=4)
c      close(unit=5)
c
c
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

      nouter = 2**10
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
        do k = 1,ninner*nbodies
          rhs(k+2*nouter) = 0.d0
          rhs(k+2*nouter+ninner) = 0.d0
        enddo
      endif
c     constant background flow on outer walls on no slip on obstacle

      do k = 1,nbodies
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+1) = 0.d0
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+2) = 0.d0
        rhs(2*nouter+2*ninner*nbodies+3*(k-1)+3) = 0.d0
      enddo
c     extra terms need for rotlet and stokeslet conditions

      end

c***********************************************************************
      subroutine solveBIE(ninner,nbodies,nouter,den,rhs,
     $    gmwork,lrwork,igwork,liwork,maxl,ifmm)
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

c     DEBUG
c      twopi = 8.d0*datan(1.d0)
c      dtheta = twopi/dble(ninner)
c      do k = 1,ninner
c        theta = dble(k-1)*dtheta
c        rhs(2*nouter + k) = dexp(dcos(theta))
c        rhs(2*nouter + k + ninner) = dexp(dsin(theta))
c      enddo

      do k = 1,2*nouter+2*ninner*nbodies+3*nbodies
        den(k) = rhs(k)
      enddo
c     initial guess

c      DEBUG
c      ntotal = 2*nouter + 2*ninner*nbodies + 3*nbodies 
c      call matvec_DLP(ntotal,den,vel1,nelt,ia,ja,a,isym)
c      call matvec_DLP_fmm(ntotal,den,vel2,nelt,ia,ja,a,isym)
c      do k = 1,ntotal
cc        if (abs(vel1(k) - vel2(k)) .ge. 1.d-10) then
c          print*,k,vel1(k)-vel2(k)
cc        endif
c      enddo


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
c     use GMRES to find the density function of the double layer
c     potential



      end


c***********************************************************************
      subroutine matvec_DLP(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**15)
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

      return
      end
c***********************************************************************
      subroutine matvec_DLP_fmm(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**15)
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

      parameter (nmax = 2**15)

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

c       START OF TARGET POINTS == OBSTCLE isou
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
      subroutine compute_shear_stress(ninner,nbodies,px,py,
     $    E11,E12,E22,shear_stress)
      implicit real*8 (a-h,o-z)

      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension E11(ninner*nbodies)
      dimension E12(ninner*nbodies)
      dimension E22(ninner*nbodies)
      dimension shear_stress(ninner*nbodies)

      dimension tractionx(ninner),tractiony(ninner)


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

c      open(unit=1,file='output/shear_stress.dat')
c      write(1,1000) shear_stress
c      close(1)

c 1000 format(E25.16)


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
      subroutine eval_velocity(ninner,nbodies,x,y,centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    ntargets,xtar,ytar,xxtar,yytar,utar,vtar,press_tar)
c     Compute the velocity and pressure at a set of target points xtar
c     and ytar.  Target points must be sufficiently far away from the
c     each boundary since no near-singular integration is used.  It is
c     just the vanilla trapezoid rule
      implicit real*8 (a-h,o-z)

      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension px(ninner*nbodies),py(ninner*nbodies)
      dimension speed(ninner*nbodies)
      dimension centerx(nbodies),centery(nbodies)
      dimension xouter(nouter),youter(nouter)
      dimension px0(nouter),py0(nouter)
      dimension speed0(nouter)
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
      dimension denx(max(ninner,nouter)),deny(max(ninner,nouter))

      dimension xxtar(ntargets),yytar(ntargets)
      dimension xtar(*),ytar(*)
      dimension utar(*),vtar(*)
      dimension press_tar(*)

      do j = 1,ntargets
        xtar(j) = xxtar(j)
        ytar(j) = yytar(j)
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




