      program stokesTracer
      implicit real*8 (a-h,o-z)
      external  f
      parameter (ninner = 2**8)
      parameter (nb = 20)
      parameter (nouter = 2**10)
      parameter (ntracer = 1000)
c      parameter (maxbodies = 10)

c      parameter (ninnc = 2**8,nbeta = 1,ninner = nbeta*ninnc)      
c      parameter (noutc = 2**8, nouter = nbeta*noutc)      
      
c      dimension centerx(maxbodies),centery(maxbodies)
c      dimension radius(maxbodies),phi(maxbodies)
      dimension x(ninner*nb),y(ninner*nb)
      dimension den(2*ninner*nb + 3*nb + 2*nouter)
c      dimension shear_stress(ninner*nbodies)
c      dimension pressure(ninner*nbodies)
c      dimension drag(2*nbodies)

      dimension xtrac(ntracer)
      dimension ytrac(ntracer)
      dimension xold(ntracer)
      dimension yold(ntracer)
      dimension x1(ntracer)      
      dimension y1(ntracer)
      dimension x2(ntracer)      
      dimension y2(ntracer)
      dimension x3(ntracer)      
      dimension y3(ntracer)      
      dimension utrac(ntracer)
      dimension vtrac(ntracer)
      dimension u1(ntracer)
      dimension v1(ntracer)
      dimension u2(ntracer)      
      dimension v2(ntracer)  
      dimension u3(ntracer)      
      dimension v3(ntracer)      
      dimension press_trac(ntracer)
      dimension vort_trac(ntracer)
      
      dimension uin(ntracer)
      dimension vin(ntracer)
      dimension xin(ntracer)
      dimension yin(ntracer)


      pi = 4.d0*datan(1.d0)
      twopi = 2.d0*pi
      dtheta = twopi/dble(ninner)
      ibary=1
      
c     the initial tracers      
      xmin = -1.d0
      xmax = 1.d0
      ymin = -0.95d0
      ymax = 0.95d0
      dt = 0.001d0
c      nx = ntargets
      ny = ntracer
      dx = (xmax - xmin)/dble(nx-1)
      dy = (ymax - ymin)/dble(ny-1)

      icount = 0
        do k = 1,ny 
          icount = icount + 1 
          xtrac(icount) = xmin 
          ytrac(icount) = ymin + dble(k-1)*dy
        enddo
        
        
      
c     input: the density function on boundary of bodies and wall

c      open(unit=1,file='output/den.dat')
c      open(unit=2,file='output/x.dat')
c      open(unit=3,file='output/y.dat')
c      
c      do k = 1,2*ninner*nbodies + 2*nouter + 3*nbodies      
c        read(1,1000) den(k)      
c      enddo
c      
c      do k = 1,ninner*nbodies
c        read(2,1000) x(k)
c        read(3,1000) y(k)
c      enddo      
c 
c      close(unit=1) 
c      close(unit=2) 
c      close(unit=3)
c

c      
!        open(unit=1,file='den/den100.txt')
!        open(unit=2,file='geo/geom100.txt')
!
!        read(2,2000) nbd
!        print *, 'n=',nbd
!                
!        do k = 1,2*ninner*nbodies + 2*nouter + 3*nbodies      
!          read(1,1000) den(k)      
!        enddo
!        do j = 1,nbodies
!          do k = 1,ninner
!            read(2,1000) x(k+(j-1)*ninner)
!          enddo
!          do k = 1,ninner
!            read(2,1000) y(k+(j-1)*ninner)
!          enddo
!        enddo      
!        close(unit=1) 
!        close(unit=2)

c     the input of gemo data        
      open(unit=1,file='geo_20b_dense/geom0001.dat')
c      open(unit=1,file='geo_100b/geom0100.dat')      
      do i=1,3
        read(1,*)
      enddo
      read(1,1020) nbodies      
      if ( nbodies .gt. 0) then
        do j=1,nbodies
          do i=1,256+3
            read(1,*)
          enddo      
          do k = 1,ninner
            read(1,1000) x(k+(j-1)*ninner)
          enddo
          do k = 1,ninner
            read(1,1000) y(k+(j-1)*ninner)
          enddo
        enddo
      endif
      close(unit=1)
      
      
c     the input of density data   
      open(unit=2,file='den_20b_dense/density0001.dat')
c      open(unit=2,file='den_100b/density0100.dat')
      do k=1, 2*ninner*nbodies + 3*nbodies + 2*nouter
        read(2,*) den(k)
      enddo
      close(unit=2)   
      
      
        

c     the number of iteration of the tracers 
      iter=5000
      irec=0
c     using the last tracer position from previous iteration as the new initial tracer
c     if the previos iteration is not enough
      if (irec . eq. 1) then
        open(unit=4,file='1/xtracer_final1.txt')    
        open(unit=5,file='1/ytracer_final1.txt')
        do k = 1,ntracer
          read(4,1000) xtrac(k)     
          read(5,1000) ytrac(k)
        enddo     
        close(unit=4) 
        close(unit=5)
      else
        open(unit=4,file='1/xtracer1.txt',
     $   access='append')  
        open(unit=5,file='1/ytracer1.txt',
     $   access='append')
        do k = 1,ntracer
          write(4,1000) xtrac(k)     
          write(5,1000) ytrac(k)
        enddo           
        close(unit=4)      
        close(unit=5)     
      endif
      
c     using the density funtion to calculate the velocities of the tracers 
c     at each iteration   
      do i = 1, iter
c        print *, xtrac(1)
        call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
     $    ntracer,xtrac,ytrac,utrac,vtrac,press_trac,vort_trac)
     
      open(unit=4,file='1/xtracer1.txt',access='append')
c     $       ,status='replace')     
      open(unit=5,file='1/ytracer1.txt',access='append')
c     $       ,status='replace')
      open(unit=6,file='1/utracer1.txt',access='append')
c     $       ,status='replace')      
      open(unit=7,file='1/vtracer1.txt',access='append')
c     $       ,status='replace')
      
      
        do k = 1,ntracer
          write(4,1000) xtrac(k)     
          write(5,1000) ytrac(k)
          write(6,1000) utrac(k)
          write(7,1000) vtrac(k)
        enddo           
     
        close(unit=4)      
        close(unit=5)      
        close(unit=6)
        close(unit=7)

c      update the tracer position by using Runge Kutta 4th method with 
c      the velocities and time step dt
        
       do k = 1,ntracer
           xold(k) = xtrac(k)
           yold(k) = ytrac(k)
c           xtrac(k)=xtrac(k)+utrac(k)*dt
c           ytrac(k)=ytrac(k)+vtrac(k)*dt
c           call rk4(xold, yold, dt, u, v, xtmp, ytmp) 
!     $           xtrac(k), ytrac(k))
           x1(k) = xold(k) + dt*utrac(k)/2.D0
           y1(k) = yold(k) + dt*vtrac(k)/2.D0
       enddo   
           
           call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
     $         ntracer,x1,y1,u1,v1,press_trac,vort_trac)
           
       do k = 1,ntracer
           x2(k) = xold(k) + dt*u1(k)/2.D0
           y2(k) = yold(k) + dt*v1(k)/2.D0
       enddo
       
           call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
     $         ntracer,x2,y2,u2,v2,press_trac,vort_trac)
           
       do k = 1,ntracer
           x3(k) = xold(k) + dt*u2(k)
           y3(k) = yold(k) + dt*v2(k)
       enddo
       
           call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
     $         ntracer,x3,y3,u3,v3,press_trac,vort_trac)
     
       do k = 1,ntracer     
           xtrac(k) = xold(k) + dt*( utrac(k) + 2.D0*u1(k) + 
     $      2.D0*u2(k) + u3(k) )/6.D0           
           ytrac(k) = yold(k) + dt*( vtrac(k) + 2.D0*v1(k) + 
     $      2.D0*v2(k) + v3(k) )/6.D0
           
           if(xtrac(k) .GT. xmax .OR. dabs(ytrac(k)) .GT. 1.d0) then
            open(unit=4,file='1/xtracer1.txt')
            open(unit=5,file='1/ytracer1.txt')
            open(unit=6,file='1/utracer1.txt')  
            open(unit=7,file='1/vtracer1.txt')
            do j = 1,ntracer
              read(4,1000) xin(j)
              read(5,1000) yin(j)            
              read(6,1000) uin(j)
              read(7,1000) vin(j)
            enddo
            close(unit=4)
            close(unit=5)
            close(unit=6)      
            close(unit=7)
            d=1000.d0
            ind=0
            do j=1,ntracer
               tmp=(uin(j)-utrac(k))**2.d0+(vin(j)-vtrac(k))**2.d0
               if(d .gt. tmp) then
                  d=tmp
                  ind=j
                endif
            enddo
             xtrac(k) = xin(ind)
             ytrac(k) = yin(ind)
           endif
       enddo
      enddo     

c     record the last position of tracer if the iteration is not enough        
        open(unit=8,file='1/xtracer_final1.txt')
c     $       ,status='replace')     
        open(unit=9,file='1/ytracer_final1.txt')
c     $       ,status='replace')       
        do k = 1,ntracer
          write(8,1000) xtrac(k)     
          write(9,1000) ytrac(k)          
        enddo
        close(unit=8)
        close(unit=9)          
        
     
     
 1000 format(E25.16)     
 1010 format(ES20.4)
 1020 format(I4)
 2000 format(I3) 
 
 
      end 
      
!c*********************************************************************      
!      subroutine rk4 ( x0, y0, dt, u0, v0, xtr, ytr)
!
!c*********************************************************************
!c
!cc RK4 takes one Runge-Kutta step.
!c
!c  Discussion:
!c
!c    It is assumed that an initial value problem, of the form
!c
!c      du/dt = f ( t, u )
!c      u(t0) = u0
!c
!c    is being solved.
!c
!c    If the user can supply current values of t, u, a stepsize dt, and a
!c    function to evaluate the derivative, this function can compute the
!c    fourth-order Runge Kutta estimate to the solution at time t+dt.
!c
!c  Licensing:
!c
!c    This code is distributed under the GNU LGPL license. 
!c
!c  Modified:
!c
!c    09 October 2013
!c
!c  Author:
!c
!c    John Burkardt
!c
!c  Parameters:
!c
!c    Input, double precision T0, the current time.
!c
!c    Input, double precision U0, the solution estimate at the current time.
!c
!c    Input, double precision DT, the time step.
!c
!c    Input, external F, a subroutine of the form 
!c      subroutine f ( t, u, uprime ) 
!c    which evaluates the derivative uprime given the time T and
!c    solution vector U.
!c
!c    Output, double precision U, the fourth-order Runge-Kutta solution 
!c    estimate at time T0+DT.
!c
!      implicit real*8 (a-h,o-z)      
!      common /geometry/x,y,ninner,nbodies,den
!      common /wall/ nouter
!c      common /tracer/ ntracer
!      common ibary
!
!      double precision dt
!      external f
!      double precision u0
!      double precision u1
!      double precision u2
!      double precision u3
!      double precision v0
!      double precision v1
!      double precision v2
!      double precision v3      
!c      double precision t0
!c      double precision t1
!c      double precision t2
!c      double precision t3
!      double precision xtr
!      double precision x0
!      double precision x1
!      double precision x2
!      double precision x3
!      double precision ytr
!      double precision y0
!      double precision y1
!      double precision y2
!      double precision y3      
!c
!c  Get four sample values of the derivative.
!c
!
!
!c      t1 = t0 + dt / 2.0D+00
!      x1 = x0 + dt*u0/2.0D+00
!      y1 = y0 + dt*v0/2.0D+00
!      
!      call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
!     $    1,x1,y1,u1,v1,press_trac,vort_trac)
!
!c      t2 = t0 + dt / 2.0D+00
!      x2 = x0 + dt*u1/2.0D+00
!      y2 = y0 + dt*v1/2.0D+00
!      
!      call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
!     $    1,x2,y2,u2,v2,press_trac,vort_trac)
!
!c      t3 = t0 + dt
!      x3 = x0 + dt*u2
!      y3 = y0 + dt*v2
!      
!      call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
!     $    1,x3,y3,u3,v3,press_trac,vort_trac)
!c
!c  Combine them to estimate the solution U at time T1.
!c
!      xtr = x0 + dt*( u0 + 2.0D+00*u1 + 2.0D+00*u2 + u3 )/6.0D+00
!      ytr = y0 + dt*( v0 + 2.0D+00*v1 + 2.0D+00*v2 + v3 )/6.0D+00
!
!      return
!      end
!      
cc****************************************************************      
c      subroutine f (xtrac, ytrac, u, v )     
cc****************************************************************
c      implicit real*8 (a-h,o-z)
c      
c      
c      common /geometry/x,y,ninner,nbodies,den
c      common /wall/ nouter
c      common /tracer/ ntracer
c      common ibary
c      
c      dimension xtrac(ntracer)
c      dimension ytrac(ntracer)
c      dimension u(ntracer)
c      dimension v(ntracer)
c
c      call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
c     $    ntracer,xtrac,ytrac,u,v,press_trac,vort_trac)
c
c      return
c      end
