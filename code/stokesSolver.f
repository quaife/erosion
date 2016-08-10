      subroutine stokesSolver(nninner,xx,yy,shear_stress)
c     Input x and y coordinates and return the shear_stress on the
c     boundary.  Outer wall is used to enclose the inner boundary so
c     that Stokes paradox is avoided
      implicit real*8 (a-h,o-z)

      dimension xx(nninner),yy(nninner)
c     x and y coordinates of obstacle

      parameter (nmax = 2**15)
c     max points on the boundary of the obstacle      
      parameter (maxl = 400, liwork = 30)
      parameter (lrwork = 10 + nmax*(maxl+6) + 
     $     maxl*(maxl+3))
c     maximum size of workspaces for GMRES

      parameter(nx = 20,ny = 20)
c     number of points in x and y direction where velocity will be
c     computed

      dimension x(nmax),y(nmax)
      dimension px(nmax), py(nmax)
c     x and y coordinates of the normal of the obstacle
      dimension cur(nmax), speed(nmax)
c     Jacobian and curvature of the geometry

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
      dimension E11(2*nninner),E12(2*nninner),E22(2*nninner)
c     components of the deformation gradient on the obstacle
      dimension shear_stress(nninner)
c     shear stress on the obstacle

      dimension xtar(nx,ny),ytar(nx,ny)
      dimension utar(nx,ny),vtar(nx,ny)

      common /geometry/ x,y,centerx,centery,px,py,cur,speed,ninner
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter
c     global variables for the geometry which are needed by the external
c     matvec routine in gmres

      ninner = nninner
c     can't declare input variables into a common field, so need new
c     variable name for x,y,ninner
      call inner_geometry(ninner,x,y,xx,yy,px,py,cur,speed)

      call outer_geometry(nouter,xouter,youter,px0,py0,cur0,speed0)
c     load geometry of initial shape

      open(unit=1,file='output/xinner.dat')
      open(unit=2,file='output/yinner.dat')
      open(unit=3,file='output/xouter.dat')
      open(unit=4,file='output/youter.dat')
      do k = 1,ninner
        write(1,1000) xx(k)
        write(2,1000) yy(k)
      enddo
      do k = 1,nouter
        write(3,1000) xouter(k)
        write(4,1000) youter(k)
      enddo
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)

      call bd_condition(ninner,x,y,nouter,xouter,youter,rhs)
c     load boundary condition

      call solveBIE(ninner,nouter,den,rhs,
     $      gmwork,lrwork,igwork,liwork,maxl)
c     solve for the density function with GMRES

      call deformation_on_boundary(ninner,x,y,centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor on the boundary of the interface

      call compute_shear_stress(ninner,px,py,E11,E12,E22,
     $    shear_stress)
c     Use the deformation tensor to compute the shear stress


      call eval_velocity(ninner,x,y,centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    nx,ny,xtar,ytar,utar,vtar)

 1000 format(E25.16)


      end


c***********************************************************************
      subroutine inner_geometry(ninner,x,y,xx,yy,px,py,cur,speed)
      implicit real*8 (a-h,o-z)

      dimension xx(ninner),yy(ninner)
      dimension x(*),y(*)
      dimension px(*), py(*)
      dimension cur(*), speed(*)


      do k = 1,ninner
        x(k) = xx(k)
        y(k) = yy(k)
      enddo

      call filter_and_derivs(ninner,x,y,px,py,cur,speed)
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
        x(k) = r*dcos(theta)
        y(k) = r*dsin(theta)
      enddo
c     outer geometry

      call filter_and_derivs(nouter,x,y,px,py,cur,speed)
      do k = 1,nouter
        px(k) = -1.d0*px(k)
        py(k) = -1.d0*py(k)
        cur(k) = -1.d0*cur(k)
      enddo
c     normal and curvatures need to point outwards



      end
        

c***********************************************************************
      subroutine bd_condition(ninner,x,y,nouter,xouter,youter,rhs)
c     Boundary condition coming from the far field flow
      implicit real*8 (a-h,o-z)

      dimension x(ninner),y(ninner)
      dimension xouter(nouter),youter(nouter)
      dimension den(2*nouter + 2*ninner+3)
      dimension rhs(2*nouter + 2*ninner+3)

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

        do k=1,ninner
          rx = x(k) - cx
          ry = y(k) - cy
          rdots = rx*sto1 + ry*sto2
          rho2 = rx**2.d0 + ry**2.d0

          rhs(k+2*nouter) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx) +
     $      rot*ry/rho2
          rhs(k+2*nouter+ninner) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry) -
     $      rot*rx/rho2
        enddo
      endif
c     velocity coming from the r \otimes r part of the Stokes flow

      if (0 .eq. 1) then
        do k = 1,nouter
          call bgFlow(1,xouter(k),youter(k),rhs(k),rhs(k+nouter))
        enddo
        do k = 1,ninner
          call bgFlow(1,x(k),y(k),
     $        rhs(k+2*nouter),rhs(k+2*nouter+ninner))
        enddo
      endif
c     generic boundary condition

      if (1 .eq. 1) then
        do k = 1,nouter
          call bgFlow(1,xouter(k),youter(k),rhs(k),rhs(k+nouter))
        enddo
        do k = 1,ninner
          rhs(k+2*nouter) = 0.d0
          rhs(k+2*nouter+ninner) = 0.d0
        enddo
      endif
c     constant flow u=[1,0] on outer walls on no slip on obstacle

      rhs(2*nouter+2*ninner+1) = 0.d0
      rhs(2*nouter+2*ninner+2) = 0.d0
      rhs(2*nouter+2*ninner+3) = 0.d0
c     extra terms need for rotlet and stokeslet conditions

      end

c***********************************************************************
      subroutine solveBIE(ninner,nouter,den,rhs,
     $    rwork,lrwork,iwork,liwork,maxl)
c     Solve the boundary integral equation with GMRES
      implicit real*8 (a-h,o-z)

      external msolve_DLP,matvec_DLP

      dimension den(2*nouter+2*ninner+3),rhs(2*nouter+2*ninner+3)
c     leave room for stokeslets and rotlets at the end of den

      dimension rwork(lrwork),iwork(liwork)
c     gmres workspaces

      itol = 0
      tol = 1.d-8
      isym = 0
      iwork(1) = maxl
      do i=2,liwork
        iwork(i) = 0
      enddo
c     paramters for DGMRES

      iwork(4) = 0
c     preconditioner flag

      iwork(5) = -1
c     restart flag

      do k = 1,2*nouter+2*ninner+3
        den(k) = rhs(k)
      enddo
c     initial guess

      call DGMRES(2*nouter+2*ninner+3,rhs,den,nelt,ia,ja,a,isym,
     $    matvec_DLP,msolve_DLP,itol,tol,itmax,iter,err,
     $    ierr,6,sb,sx,rwork,lrwork,iwork,liwork,rw,iw)
c     use GMRES to find the density function of the double layer
c     potential

      end


c***********************************************************************
      subroutine matvec_DLP(ntotal,den,vel,nelt,ia,ja,a,isym)
c     matrix vector multiplication routine for the double-layer
c     potential
      implicit real*8 (a-h,o-z)

      parameter (nmax = 2**15)

      dimension den(ntotal)
      dimension vel(ntotal)

      common /geometry/ x,y,centerx,centery,px,py,cur,speed,ninner
      common /wall/ xouter,youter,px0,py0,cur0,speed0,nouter

      dimension x(nmax),y(nmax)
      dimension px(nmax),py(nmax)
      dimension cur(nmax),speed(nmax)

      dimension xouter(nmax),youter(nmax)
      dimension px0(nmax),py0(nmax)
      dimension cur0(nmax),speed0(nmax)

      dimension denx(nmax),deny(nmax)
      dimension ux(nmax),uy(nmax)

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi

      do k = 1,2*ninner + 2*nouter + 3
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
      do k=1,nouter
        ux(k) = 0.d0
        uy(k) = 0.d0

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

c     START OF TARGET POINTS == OBSTACLE
      do k=1,ninner
        ux(k) = 0.d0
        uy(k) = 0.d0

        do j=1,nouter
          rx = x(k) - xouter(j)
          ry = y(k) - youter(j)
          rdotn = rx*px0(j) + ry*py0(j)
          rdotden = rx*denx(j) + ry*deny(j)
          rho2 = rx**2.d0 + ry**2.d0

          ux(k) = ux(k) + 
     $        rdotn*rdotden/rho2/rho2*rx*speed0(j)
          uy(k) = uy(k) + 
     $        rdotn*rdotden/rho2/rho2*ry*speed0(j)
        enddo
        ux(k) = ux(k)*twopi/dble(nouter)/pi
        uy(k) = uy(k)*twopi/dble(nouter)/pi
      enddo

      do k=1,ninner
        vel(k + 2*nouter) = vel(k + 2*nouter) + ux(k)
        vel(k + 2*nouter + ninner) = vel(k + 2*nouter + ninner) + uy(k)
      enddo
c     END OF TARGET POINTS == OBSTACLE

      do k = 1,nouter
        vel(k) = vel(k) - 5.d-1*denx(k)
        vel(k + nouter) = vel(k + nouter) - 5.d-1*deny(k)
      enddo
c     add in the jump condition
c     END OF SOURCE POINTS == OUTER WALL

c************************************************************

c     START OF SOURCE POINTS == OBSTACLE
      do k = 1,ninner
        denx(k) = den(k + 2*nouter)
        deny(k) = den(k + 2*nouter + ninner)
      enddo
c     density function due to obstacle

c     START OF TARGET POINTS == OBSTACLE
      do k=1,ninner
        ux(k) = 0.d0
        uy(k) = 0.d0

        do j=1,ninner
          if (j .ne. k) then
            rx = x(k) - x(j)
            ry = y(k) - y(j)
            rdotn = rx*px(j) + ry*py(j)
            rdotden = rx*denx(j) + ry*deny(j)
            rho2 = rx**2.d0 + ry**2.d0

            ux(k) = ux(k) + 
     $          rdotn*rdotden/rho2/rho2*rx*speed(j)
            uy(k) = uy(k) + 
     $          rdotn*rdotden/rho2/rho2*ry*speed(j)
          endif
        enddo
        ux(k) = ux(k)*twopi/dble(ninner)/pi
        uy(k) = uy(k)*twopi/dble(ninner)/pi
      enddo

      do k=1,ninner
        tdotden = py(k)*denx(k) - px(k)*deny(k)
        vel(k + 2*nouter) = vel(k + 2*nouter) +
     $        ux(k) - cur(k)*tdotden*py(k)*
     $        speed(k)/dble(ninner)
        vel(k + 2*nouter + ninner) = vel(k + 2*nouter + ninner) + 
     $        uy(k) + cur(k)*tdotden*px(k)*
     $        speed(k)/dble(ninner)
      enddo
c     Add in correction at diagonal term that involves the curvature
c     END OF TARGET POINTS == OBSTACLE

c     START OF TARGET POINTS == OUTER WALL
      do k=1,nouter
        ux(k) = 0.d0
        uy(k) = 0.d0

        do j=1,ninner
          rx = xouter(k) - x(j)
          ry = youter(k) - y(j)
          rdotn = rx*px(j) + ry*py(j)
          rdotden = rx*denx(j) + ry*deny(j)
          rho2 = rx**2.d0 + ry**2.d0

          ux(k) = ux(k) + 
     $        rdotn*rdotden/rho2/rho2*rx*speed(j)
          uy(k) = uy(k) + 
     $        rdotn*rdotden/rho2/rho2*ry*speed(j)
        enddo
        ux(k) = ux(k)*twopi/dble(ninner)/pi
        uy(k) = uy(k)*twopi/dble(ninner)/pi
      enddo

      do k=1,nouter
        vel(k) = vel(k) + ux(k)
        vel(k + nouter) = vel(k + nouter) + uy(k)
      enddo
c     END OF TARGET POINTS == OUTER WALL

      do k = 1,ninner
        vel(k + 2*nouter) = vel(k + 2*nouter) - 5.d-1*denx(k)
        vel(k + 2*nouter + ninner) = vel(k + 2*nouter + ninner) - 
     $      5.d-1*deny(k)
      enddo
c     add in the jump condition
c     END OF TARGET POINTS == OUTER WALL

c************************************************************

c     START OF SOURCE POINTS == ROTLETS AND STOKESLETS
      sto1 = den(2*nouter+2*ninner+1)
      sto2 = den(2*nouter+2*ninner+2)
      rot = den(2*nouter+2*ninner+3)

c     START OF TARGET POINTS == OUTER WALL
      do k = 1,nouter
        ux(k) = 0.d0
        uy(k) = 0.d0
        rx = xouter(k) - centerx
        ry = youter(k) - centery
        rho2 = rx**2.d0 + ry**2.d0
        rdots = rx*sto1 + ry*sto2
        ux(k) = 5.d-1/twopi*
     $      (-5.d-1*log(rho2)*sto1 + rdots/rho2*rx)
        uy(k) = 5.d-1/twopi*
     $      (-5.d-1*log(rho2)*sto2 + rdots/rho2*ry)
        ux(k) = ux(k) + rot*ry/rho2
        uy(k) = uy(k) - rot*rx/rho2
      enddo
c     stokeslets and rotlets contribute to the velocity

      do k=1,nouter
        vel(k) = vel(k) + ux(k)
        vel(k+nouter) = vel(k+nouter) + uy(k)
      enddo
c     END OF TARGET POINTS == OUTER WALL

c     START OF TARGET POINTS == OBSTACLE
      do k = 1,ninner
        ux(k) = 0.d0
        uy(k) = 0.d0
        rx = x(k) - centerx
        ry = y(k) - centery
        rho2 = rx**2.d0 + ry**2.d0
        rdots = rx*sto1 + ry*sto2
        ux(k) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
        uy(k) = 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
        ux(k) = ux(k) + rot*ry/rho2
        uy(k) = uy(k) - rot*rx/rho2
      enddo
c     stokeslets and rotlets contribute to the velocity

      do k=1,ninner
        vel(k + 2*nouter) = vel(k + 2*nouter) + ux(k)
        vel(k + 2*nouter + ninner) = vel(k + 2*nouter + ninner) + 
     $      uy(k)
      enddo
c     END OF TARGET POINTS == OBSTACLE

c     END OF SOURCE POINTS == ROTLETS AND STOKESLETS

c************************************************************

c     START OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c     STOKESLETS
      do k = 1,ninner
        denx(k) = den(k + 2*nouter)
        deny(k) = den(k + 2*nouter + ninner)
      enddo
c     density function due to obstacle
      vel(2*nouter+2*ninner+1) = 0.d0
      vel(2*nouter+2*ninner+2) = 0.d0
      vel(2*nouter+2*ninner+3) = 0.d0
      do k = 1,ninner
        vel(2*nouter+2*ninner+1) = vel(2*nouter+2*ninner+1) - 
     $      denx(k)*speed(k)
        vel(2*nouter+2*ninner+2) = vel(2*nouter+2*ninner+2) - 
     $      deny(k)*speed(k)
        vel(2*nouter+2*ninner+3) = vel(2*nouter+2*ninner+3) -
     $        (denx(k)*y(k)-deny(k)*x(k))*speed(k)
      enddo
      vel(2*nouter+2*ninner+1) = vel(2*nouter+2*ninner+1)/twopi * 
     $    (twopi/dble(ninner))
      vel(2*nouter+2*ninner+2) = vel(2*nouter+2*ninner+2)/twopi * 
     $    (twopi/dble(ninner))
      vel(2*nouter+2*ninner+3) = vel(2*nouter+2*ninner+3)/twopi * 
     $    (twopi/dble(ninner))
      vel(2*nouter+2*ninner+1) = vel(2*nouter+2*ninner+1) + sto1
      vel(2*nouter+2*ninner+2) = vel(2*nouter+2*ninner+2) + sto2
      vel(2*nouter+2*ninner+3) = vel(2*nouter+2*ninner+3) + rot
c     END OF INTEGRALS OF DENSITY FUNCTION BEING EQUAL TO ROTLETS AND
c     STOKESLETS

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

      dimension r(nn),z(nn)

      do i = 1,nn
        r(i) = z(i)
      enddo
c     no preconditioner for now

      return
      end


c***********************************************************************
      subroutine deformation_on_boundary(ninner,x,y,centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    E11,E12,E22)
c     compute the deformation tensor at points on the boundary of the
c     obstacle
      implicit real*8 (a-h,o-z)

      dimension x(ninner),y(ninner)
      dimension px(ninner),py(ninner)
      dimension speed(ninner)
      dimension xouter(nouter),youter(nouter)
      dimension px0(nouter),py0(nouter)
      dimension speed0(nouter)
      dimension den(2*nouter + 2*ninner + 3)

      dimension E11(ninner),E12(ninner),E22(ninner)

      complex *16 zden(ninner)
      real*8 wsave(4*ninner+15)
      complex *16 eye
      real *8 denx(max(ninner,nouter)),deny(max(ninner,nouter))

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      eye = (0.d0,1.d0)

c     Contribution from the density function on the solid wall
      do j=1,nouter
        denx(j) = den(j)
        deny(j) = den(j+nouter)
      enddo
c     Density function defined on the outer geometry

      do j=1,ninner
        E11(j) = 0.d0
        E12(j) = 0.d0
        E22(j) = 0.d0
c       Initialize contributions to be 0
        do k=1,nouter
          rx = x(j) - xouter(k)
          ry = y(j) - youter(k)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px0(k) + ry*py0(k)
          routn = rx*py0(k) + ry*px0(k)
          rdotden = rx*denx(k) + ry*deny(k)
          routden = rx*deny(k) + ry*denx(k)

          E11(j) = E11(j) + 5.d-1*
     $        (2.d0*rdotn*rdotden/rho2/rho2 + 
     $        2.d0*rdotden/rho2/rho2*rx*px0(k) +
     $        2.d0*rdotn/rho2/rho2*rx*denx(k) -
     $        8.d0*rdotn*rdotden/rho2**3.d0*rx*rx)*
     $        speed0(k)*twopi/dble(nouter)/pi

          E12(j) = E12(j) + 5.d-1*
     $        (rdotden*routn/rho2/rho2 +
     $        rdotn*routden/rho2/rho2 -
     $        8.d0*rdotn*rdotden*rx*ry/rho2**3.d0)*
     $        speed0(k)*twopi/dble(nouter)/pi

          E22(j) = E22(j) + 5.d-1*
     $        (2.d0*rdotn*rdotden/rho2/rho2 + 
     $        2.d0*rdotden/rho2/rho2*ry*py0(k) +
     $        2.d0*rdotn/rho2/rho2*ry*deny(k) -
     $        8.d0*rdotn*rdotden/rho2**3.d0*ry*ry)*
     $        speed0(k)*twopi/dble(nouter)/pi
        enddo
      enddo

      do j=1,ninner
        denx(j) = den(j+2*nouter)
        deny(j) = den(j+2*nouter+ninner)
      enddo
c     Density function defined on the inner geometry

c     Contribution from the density function on the obstacle itself
      do j = 1,ninner,2
        do k = 2,ninner,2
          rx = x(k) - x(j)
          ry = y(k) - y(j)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px(k) + ry*py(k)
          routn = rx*py(k) + ry*px(k)

          sx = denx(k) - denx(j)
          sy = deny(k) - deny(j)
c         subtract of the density function evaluated at the target point
c         to reduce the singularity to 1/r.  Then, odd-even integration
c         converges to the PV integral
          rdots = rx*sx + ry*sy
          routs = rx*sy + ry*sx

          E11(j) = E11(j) + 2.d0*5.d-1*
     $        (2.d0*rdotn*rdots/rho2/rho2 + 
     $        2.d0*rdots/rho2/rho2*rx*px(k) +
     $        2.d0*rdotn/rho2/rho2*rx*sx -
     $        8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
     $        speed(k)*twopi/dble(ninner)/pi

          E12(j) = E12(j) + 2.d0*5.d-1*
     $        (rdots*routn/rho2/rho2 +
     $        rdotn*routs/rho2/rho2 -
     $        8.d0*rdotn*rdots*rx*ry/rho2**3.d0)*
     $        speed(k)*twopi/dble(ninner)/pi

          E22(j) = E22(j) + 2.d0*5.d-1*
     $        (2.d0*rdotn*rdots/rho2/rho2 + 
     $        2.d0*rdots/rho2/rho2*ry*py(k) +
     $        2.d0*rdotn/rho2/rho2*ry*sy -
     $        8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
     $        speed(k)*twopi/dble(ninner)/pi
        enddo
c       need to multiply by 2 since the grid spacing is twice as large
      enddo
c     compute deformation tensor on odd indexed terms

      do j = 2,ninner,2
        do k = 1,ninner,2
          rx = x(k) - x(j)
          ry = y(k) - y(j)
          rho2 = rx**2.d0 + ry**2.d0
          rdotn = rx*px(k) + ry*py(k)
          routn = rx*py(k) + ry*px(k)

          sx = denx(k) - denx(j)
          sy = deny(k) - deny(j)
c         subtract of the density function evaluated at the target point
c         to reduce the singularity to 1/r.  Then, odd-even integration
c         converges to the PV integral
          rdots = rx*sx + ry*sy
          routs = rx*sy + ry*sx

          E11(j) = E11(j) + 2.d0*5.d-1*
     $        (2.d0*rdotn*rdots/rho2/rho2 + 
     $        2.d0*rdots/rho2/rho2*rx*px(k) +
     $        2.d0*rdotn/rho2/rho2*rx*sx -
     $        8.d0*rdotn*rdots/rho2**3.d0*rx*rx)*
     $        speed(k)*twopi/dble(ninner)/pi

          E12(j) = E12(j) + 2.d0*5.d-1*
     $        (rdots*routn/rho2/rho2 +
     $        rdotn*routs/rho2/rho2 -
     $        8.0d0*rdotn*rdots*rx*ry/rho2**3.d0)*
     $        speed(k)*twopi/dble(ninner)/pi

          E22(j) = E22(j) + 2.d0*5.d-1*
     $        (2.d0*rdotn*rdots/rho2/rho2 + 
     $        2.d0*rdots/rho2/rho2*ry*py(k) +
     $        2.d0*rdotn/rho2/rho2*ry*sy -
     $        8.d0*rdotn*rdots/rho2**3.d0*ry*ry)*
     $        speed(k)*twopi/dble(ninner)/pi
        enddo
c       need to multiply by 2 since the grid spacing is twice as large
      enddo
c     compute deformation tensor on odd indexed terms

      call DCFFTI(ninner,wsave)
      do k=1,ninner 
        zden(k) = denx(k) + eye*deny(k)
      enddo
      call fourierDiff(ninner,zden,wsave)
c     real part of zden is parameter derivative of the first component
c     of the density function and the complex part is the derivative of
c     the second component of the density function

      do j=1,ninner
        tx = -py(j)
        ty = px(j)
        dsdtx = dreal(zden(j))/speed(j)
        dsdty = dimag(zden(j))/speed(j)
        dsdt_dot_tau = dsdtx*tx + dsdty*ty
        E11(j) = E11(j) + 5.d-1*
     $      dsdt_dot_tau*(tx**2.d0 - ty**2.d0)
        E12(j) = E12(j) + 5.d-1*
     $      dsdt_dot_tau*2.d0*tx*ty
        E22(j) = E22(j) + 5.d-1*
     $      dsdt_dot_tau*(-tx**2.d0 + ty**2.d0)
      enddo
c     Add in jump conditions

c     CONTRIBUTION FROM ROTLETS AND STOKESLETS
      sto1 = den(2*nouter+2*ninner+1)
      sto2 = den(2*nouter+2*ninner+2)
      rot = den(2*nouter+2*ninner+3)
      do j = 1,ninner
        rx = x(j) - centerx
        ry = y(j) - centery
        rho2 = rx**2.d0 + ry**2.d0
        rdots = rx*sto1 + ry*sto2
        E11(j) = E11(j) + 5.d-1/twopi*
     $      (rdots/rho2 - 2*rx*rx*rdots/rho2/rho2)
        E22(j) = E22(j) + 5.d-1/twopi*
     $      (rdots/rho2 - 2*ry*ry*rdots/rho2/rho2)
        E12(j) = E12(j) - 1.d0/twopi*
     $      (rdots/rho2/rho2*rx*ry)
c       stokeslet contribution

        E11(j) = E11(j) - 2.d0*rot*rx*ry/rho2/rho2
        E12(j) = E12(j) + rot*(rx**2.d0 - ry**2.d0)/rho2/rho2
        E22(j) = E22(j) + 2.d0*rot*rx*ry/rho2/rho2
c       rotlet contribution
      enddo
c     Add in contribution from rotlets and stokeslets

      open(unit=1,file='output/den.dat')
      write(1,1000) den
      close(1)

 1000 format(E25.16)

      end

c***********************************************************************
      subroutine compute_shear_stress(ninner,px,py,E11,E12,E22,
     $    shear_stress)
      implicit real*8 (a-h,o-z)

      dimension px(ninner),py(ninner)
      dimension E11(ninner),E12(ninner),E22(ninner)
      dimension shear_stress(ninner)

      dimension tractionx(ninner),tractiony(ninner)

      tractionx = -2.d0*E11*px - 2.d0*E12*py
      tractiony = -2.d0*E12*px - 2.d0*E22*py
c     Point normal outward as Nick uses

      shear_stress = tractionx*py - tractiony*px 
c     Finally compute shear stress

      open(unit=1,file='output/shear_stress.dat')
      write(1,1000) shear_stress
      close(1)

 1000 format(E25.16)


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
        u(k) = 1.d0
        v(k) = 0.d0
c       this flow has vanishing deformation tensor

c        u(k) = x(k)
c        v(k) = -y(k)
c       this flow does have a non-trivial deformation tensor

c        u(k) = 2.d0*y(k) - y(k)**2.d0 - x(k)
c        v(k) = y(k)
c       this flow does have a non-trivial deformation tensor

c        u(k) = x(k)/(x(k)**2.d0 + y(k)**2.d0)
c        v(k) = y(k)/(x(k)**2.d0 + y(k)**2.d0)
c       this flow does have a non-trivial deformation tensor

c        u(k) = 4.0d0-y(k)**2.d0
c        v(k) = 0.d0
c       this flow does have a non-trivial deformation tensor

c        u(k) = y(k)
c        v(k) = 0.d0
c       this flow does have a non-trivial deformation tensor
    
c        u(k) = y(k)/(x(k)**2.d0 + y(k)**2.d0)
c        v(k) = -x(k)/(x(k)**2.d0 + y(k)**2.d0)
c       single Rotlet
      enddo

      end 


c***********************************************************************
      subroutine filter_and_derivs(n,x,y,px,py,cur,speed)
      implicit real*8 (a-h,o-z)

      dimension x(n),y(n)
      dimension px(n),py(n)
      dimension cur(n),speed(n)

      real *8 Dx(n),Dy(n)
      real *8 DDx(n),DDy(n)

      complex *16 zden1(n),zden2(n)
      real*8 wsave(4*n+15)
      complex *16 eye

      eye = (0.d0,1.d0)

      zden1 = x + eye*y
      call fourierDiff(n,zden1,wsave)
      zden2 = zden1
      call fourierDiff(n,zden2,wsave)
c     real part of zden is the first derivative of the x component of
c     the position, imaginary part of zden is the first derivative of
c     the y component of the position

      Dx = dreal(zden1)
      Dy = dimag(zden1)
      speed = dsqrt(Dx**2.d0 + Dy**2.d0)
      px = -Dy/speed
      py = Dx/speed
      DDx = dreal(zden2)
      DDy = dimag(zden2)
      cur = (Dy*DDx - Dx*DDy)/speed**3.d0
        

      end



c***********************************************************************
      subroutine eval_velocity(ninner,x,y,centerx,centery,
     $    px,py,speed,nouter,xouter,youter,px0,py0,speed0,den,
     $    nx,ny,xtar,ytar,utar,vtar)
c     compute the velocity on a meshgrid
      implicit real*8 (a-h,o-z)

      dimension x(ninner),y(ninner)
      dimension px(ninner),py(ninner)
      dimension speed(ninner)
      dimension xouter(nouter),youter(nouter)
      dimension px0(nouter),py0(nouter)
      dimension speed0(nouter)
      dimension den(2*nouter + 2*ninner + 3)
      dimension xtar(nx,ny),ytar(nx,ny),utar(nx,ny),vtar(nx,ny)

c      complex *16 zden(ninner)
c      real*8 wsave(4*ninner+15)
c      complex *16 eye
      real *8 denx(max(ninner,nouter)),deny(max(ninner,nouter))

      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
c      eye = (0.d0,1.d0)

      xstart = -8.d-1
      xend = 8.d-1
      ystart = -8.d-1
      yend = 8.d-1
      dx = (xend - xstart)/dble(nx-1)
      dy = (yend - ystart)/dble(ny-1)
      do j = 1,nx
        do k = 1,ny
          xtar(j,k) = xstart + dble(j-1)*dx
          ytar(j,k) = ystart + dble(k-1)*dy
          utar(j,k) = 0.d0
          vtar(j,k) = 0.d0
        enddo
      enddo

c     Contribution from the density function on the solid wall
      do j=1,nouter
        denx(j) = den(j)
        deny(j) = den(j+nouter)
      enddo
c     Density function defined on the outer geometry

      do i = 1,nx
        do j = 1,ny
          do k = 1,nouter
            rx = xtar(i,j) - xouter(k)
            ry = ytar(i,j) - youter(k)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px0(k) + ry*py0(k)
            rdotden = rx*denx(k) + ry*deny(k)

            utar(i,j) = utar(i,j) + rdotn/rho2*rdotden/rho2*rx*
     $          speed0(k)*twopi/dble(nouter)/pi
            vtar(i,j) = vtar(i,j) + rdotn/rho2*rdotden/rho2*ry*
     $          speed0(k)*twopi/dble(nouter)/pi
          enddo
        enddo
      enddo

c     Contribution from the density function on the inner wall
      do j=1,ninner
        denx(j) = den(j+2*nouter)
        deny(j) = den(j+2*nouter+ninner)
      enddo
c     Density function defined on the inner geometry

      do i = 1,nx
        do j = 1,ny
          do k = 1,ninner
            rx = xtar(i,j) - x(k)
            ry = ytar(i,j) - y(k)
            rho2 = rx**2.d0 + ry**2.d0
            rdotn = rx*px(k) + ry*py(k)
            rdotden = rx*denx(k) + ry*deny(k)

            utar(i,j) = utar(i,j) + rdotn/rho2*rdotden/rho2*rx*
     $          speed(k)*twopi/dble(ninner)/pi
            vtar(i,j) = vtar(i,j) + rdotn/rho2*rdotden/rho2*ry*
     $          speed(k)*twopi/dble(ninner)/pi
          enddo
        enddo
      enddo

c     CONTRIBUTION FROM ROTLETS AND STOKESLETS
      sto1 = den(2*nouter+2*ninner+1)
      sto2 = den(2*nouter+2*ninner+2)
      rot = den(2*nouter+2*ninner+3)
      do i = 1,nx
        do j = 1,ny
          rx = xtar(i,j) - centerx
          ry = ytar(i,j) - centery
          rho2 = rx**2.d0 + ry**2.d0
          rdots = rx*sto1 + ry*sto2
          utar(i,j) = utar(i,j) + 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto1 + rdots/rho2*rx)
          vtar(i,j) = vtar(i,j) + 5.d-1/twopi*
     $      (-5.d-1*dlog(rho2)*sto2 + rdots/rho2*ry)
c         stokeslet contribution

          utar(i,j) = utar(i,j) + rot*ry/rho2
          vtar(i,j) = vtar(i,j) - rot*rx/rho2
c         rotlet contribution
        enddo
      enddo
c     Add in contribution from rotlets and stokeslets

      open(unit=1,file='output/xx.dat')
      open(unit=2,file='output/yy.dat')
      open(unit=3,file='output/uu.dat')
      open(unit=4,file='output/vv.dat')
      do j = 1,nx
        write(1,1000) (xtar(j,k), k=1,ny)
        write(2,1000) (ytar(j,k), k=1,ny)
        write(3,1000) (utar(j,k), k=1,ny)
        write(4,1000) (vtar(j,k), k=1,ny)
      enddo
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)

 1000 format(100(E25.16))



      end

