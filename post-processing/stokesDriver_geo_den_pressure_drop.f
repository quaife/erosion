      program stokesDriver
      implicit real*8 (a-h,o-z)

      parameter (ninner = 2**8)
      parameter (nmax = 50)
      parameter (nouter = 2**10)
      parameter (ntargets = 500)
      parameter (maxbodies = 100)

c      parameter (ninnc = 2**8,nbeta = 1,ninner = nbeta*ninnc)      
c      parameter (noutc = 2**8, nouter = nbeta*noutc)      
      
      dimension centerx(maxbodies),centery(maxbodies)
      dimension radius(maxbodies),phi(maxbodies)
      dimension x(ninner*nmax),y(ninner*nmax)
      dimension den(2*ninner*nmax + 3*nmax + 2*nouter)
c      dimension shear_stress(ninner*nbodies)
c      dimension pressure(ninner*nbodies)
c      dimension drag(2*nbodies)

!      real *8, allocatable :: x(:)
!      real *8, allocatable :: y(:)
!      real *8, allocatable :: den(:)      
!      real *8, allocatable :: shear_stress(:)
!      real *8, allocatable :: pressure(:)
!      real *8, allocatable :: drag(:)      
      
      dimension xtar(ntargets*ntargets)
      dimension ytar(ntargets*ntargets)
      dimension utar(ntargets*ntargets)
      dimension vtar(ntargets*ntargets)
      dimension press_tar(ntargets*ntargets)
      dimension vort_tar(ntargets*ntargets)
      
c     Testing for stress tensor       
c      dimension xtar_test(ninner*nbodies)
c      dimension ytar_test(ninner*nbodies)
c      dimension utar_test(ninner*nbodies)
c      dimension vtar_test(ninner*nbodies)
c      dimension press_tar_test(ninner*nbodies)
c      dimension vort_tar_test(ninner*nbodies)      
      
      
      dimension iside(ntargets*ntargets)
      dimension inear(ntargets*ntargets)
      
      character(len=1024) :: f1      
      character(len=1024) :: format_string1
      character(len=1024) :: f2      
      character(len=1024) :: format_string2
      
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


!      radius(1) = 2.d-1
!      radius(2) = 2.d-1      
!      centerx(1) = 0.5d0
!      centerx(2) = 0.5d0 - radius(1) - radius(2) - 0.001d0
!      centery(1) = 0.1d0 + 0.499d0
!      centery(2) = 0.1d0 + 0.499d0
!      if(nbodies .GT. 2) then
!        j = 2
!        radius(1) = 2.d-1
!        centerx(1) = 0.5d0
!        centery(1) = 0.1d0 + 0.499d0        
!        do n= 1,1000
!          centerx(j) = -0.799d0+rand()*0.799d0*2
!          centery(j) = -0.799d0+rand()*0.799d0*2
!          radius(j) = 2.d-1
!          do k=1,j-1
!            tmp=dsqrt((centerx(j)-centerx(k))**2+
!     &              (centery(j)-centery(k))**2)
!            print *, 'tmp=',tmp
!            if(tmp .lt. radius(1)*2.d0) go to 10
!          enddo
!          j = j + 1
!          if(j .eq. nbodies +1) EXIT
!10      enddo
!      endif

       
!!c      radius(3) = 5.d-2
!      phi(1) = 0.d0
!      phi(2) = 0.d0
!!c      phi(3) = pi/4.d0
      
c      Nphi = 128
c      dphi = twopi/Nphi
c      open(unit=21,file='output/dragx.dat')
c      open(unit=22,file='output/dragy.dat')

c      do kk = 1,Nphi
c        phi(1) = dble(kk-1)*dphi

c      do j = 1,nbodies
c        do k = 1,ninner
c          theta = dble(k-1)*dtheta
c          var_rad = radius(j)
c          x((j-1)*ninner+k) = centerx(j) + var_rad*dcos(theta)
c          y((j-1)*ninner+k) = centery(j) + var_rad*dsin(theta)
c!          x((j-1)*ninner+k) = centerx(j) + var_rad*
c!     $          (dcos(phi(j))*dcos(theta) + 
c!     $           2.d0*dsin(phi(j))*dsin(theta))
c!          y((j-1)*ninner+k) = centery(j) + var_rad*
c!     $          (-dsin(phi(j))*dcos(theta) + 
c!     $           2.d0*dcos(phi(j))*dsin(theta))
c        enddo
c      enddo


c     input
      ibary = 1
      ifmm = 1
      nbeta = 1
      maxl = 10000


c     target points
      xmin = -1.d0
      xmax = 1.d0
      ymin = -0.99d0
      ymax = 0.99d0
      nx =2
c      nx = ntargets
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
      
c      if ( ibary .eq. 1) then
c          if( ifmm .eq. 1) then 
c             open(unit=1,file='record/xtar_fmmbary.dat')          
c             open(unit=2,file='record/ytar_fmmbary.dat')          
c          else
c             open(unit=1,file='output/xtar.dat')
c             open(unit=2,file='output/ytar.dat')
c          endif
c      else
c        if( ifmm .eq. 1) then
c          open(unit=1,file='output/xtar_fmm.dat')      
c          open(unit=2,file='output/ytar_fmm.dat')      
c        else
c          open(unit=1,file='output/trap_xtar.dat')      
c          open(unit=2,file='output/trap_ytar.dat')      
c        endif
c      endif
c      
c      do k = 1,ntargets**2
c        write(1,1000) xtar(k)
c        write(2,1000) ytar(k)
c      enddo
c      close(unit=1)
c      close(unit=2)     

      format_string1 = "(I0)"
      format_string2 = "(I4.4)"
      
      iterm = 347
      do istep= 1, iterm
        print *, 'start of loop of step', istep      
       
          write (f1,format_string1) istep  
          write (f2,format_string2) istep
!        open(unit=1,file='geo_20b_dense/geom.dat')  
        open(unit=1,file='run_params347stop/geom'//trim(f2)//'.dat')
        do i=1,3
          read(1,*)
        enddo
        read(1,1020) nbodies      
!        allocate(x(ninner*nbodies))
!        allocate(y(ninner*nbodies)) 
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
        
!        allocate(den(2*ninner*nbodies + 3*nbodies + 2*nouter))
!        allocate(shear_stress(ninner*nbodies))
!        allocate(pressure(ninner*nbodies))
!        allocate(drag(2*nbodies))      

!        open(unit=2,file='den_20b_dense/density0100.dat')
        open(unit=2,file='run_params347stop/density'//trim(f2)//'.dat')
        do k=1, 2*ninner*nbodies + 3*nbodies + 2*nouter
          read(2,*) den(k)
        enddo
        close(unit=2)
  
  
!      call stokesSolver(ninner,nbodies,nouter,ifmm,ibary,maxl,
!     $       x,y,den,iter)
!c     pass in number of points and x and y coordinates and return the
!c     density function on the boundary
!      print *, 'Solver'
c      call classifyPoints(ninner,nbodies,x,y,
c     $      ntargets*ntargets,xtar,ytar,iside,inear)
c
c      call computeVorticityTargets(ninner,nbodies,nouter,
c     $      x,y,den,ntargets*ntargets,xtar,ytar,vort_tar)
c
c      call computeVelocityPressureTargets(ninner,nbodies,nouter,
c     $    x,y,den,ntargets*ntargets,xtar,ytar,utar,vtar,press_tar)
c
c  
c        if ( ibary .eq. 1 .and. ifmm .eq. 1) then
c          open(unit=1,file='record/nbody.dat',access = 'append')
c          open(unit=8,file='record/x'//trim(f1)//'.dat')
c          open(unit=9,file='record/y'//trim(f1)//'.dat')
c        elseif (ibary .eq. 1 .and. ifmm .eq. 0) then  
c          open(unit=8,file='output/x.dat')
c          open(unit=9,file='output/y.dat')
c        elseif (ibary .eq. 0 .and. ifmm .eq. 1) then    
c          open(unit=8,file='output/x_fmm.dat')      
c          open(unit=9,file='output/y_fmm.dat')
c        else   
c          open(unit=8,file='output/trap_x.dat')      
c          open(unit=9,file='output/trap_y.dat')      
c        endif
c        write(1,1020) nbodies
c        do k = 1,ninner*nbodies
c          write(8,1000) x(k)
c          write(9,1000) y(k)
c        enddo
c        close(unit=1)
c        close(unit=8)
c        close(unit=9)
  
c     output 
        ibary = 1      
!        if( ibary .eq. 1) then
!        print *, 'computeQoiTargets_Bary'
!        else 
!c      call computeQoiTargetsTrap(ninner,nbodies,nouter,x,y,den,
!c     $  ntargets*ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)
!        print *, 'computeQoiTargets_trap'      
!        endif
!        
!        
        call computeQoiTargets(ninner,nbodies,nouter,ibary,x,y,den,
     $     nx*ny,xtar,ytar,utar,vtar,press_tar,vort_tar)
       
!        call computeShearStress(ninner,nbodies,nouter,x,y,den,ibary,
!       $    shear_stress)
!c     pass in the density function and return the shear_stress
!      print *, 'computeShearStress'
!      call computePressure(ninner,nbodies,nouter,x,y,den,pressure)
!c     pass in the density function and return the pressure
!      print *, 'computePressure'
!      call computeDrag(ninner,nbodies,x,y,
!     $      shear_stress,pressure,drag)
!c     pass in the shear_stress and pressure and return the drag
!c      write(21,1000) drag(1)
!c      write(22,1000) drag(2)
!c      enddo
!c      close(unit=21)
!c      close(unit=22)
!        print *, 'computeDrag' 
!        
  
  !      if ( ibary .eq. 1 .and. ifmm .eq. 1) then
  !        open(unit=2,file='output/record/shear_stress_fmmbary.dat')
  !        open(unit=3,file='output/record/pressure_fmmbary.dat')
  !        open(unit=4,file='output/record/drag_fmmbary.dat')
  !      elseif (ibary .eq. 1 .and. ifmm .eq. 0) then  
  !        open(unit=2,file='output/shear_stress.dat')
  !        open(unit=3,file='output/pressure.dat')
  !        open(unit=4,file='output/drag.dat')
  !      elseif (ibary .eq. 0 .and. ifmm .eq. 1) then
  !        open(unit=2,file='output/shear_stress_fmm.dat') 
  !        open(unit=3,file='output/pressure_fmm.dat') 
  !        open(unit=4,file='output/drag_fmm.dat')      
  !      else 
  !        open(unit=2,file='output/trap_shear_stress.dat') 
  !        open(unit=3,file='output/trap_pressure.dat') 
  !        open(unit=4,file='output/trap_drag.dat')           
  !      endif
  !      
  !      do k = 1,ninner*nbodies
  !        write(2,1000) shear_stress(k)
  !        write(3,1000) pressure(k)
  !      enddo
  !      do k = 1,2*nbodies
  !        write(4,1000) drag(k)
  !      enddo
  !      close(unit=2)
  !      close(unit=3)
  !      close(unit=4)
  
c        if ( ibary .eq. 1) then
c          if( ifmm .eq. 1) then          
c            open(unit=3,file='record/utar'//trim(f1)//'.dat')          
c            open(unit=4,file='record/vtar'//trim(f1)//'.dat')
c            open(unit=8,file='record/press_tar'//trim(f1)//'.dat')            
            open(unit=8,
     &       file='pressure_drop/drop_press'//trim(f1)//'.dat')          
c            open(unit=9,file='record/vort_tar'//trim(f1)//'.dat')
!!            open(unit=3,file='output/utar100.dat')          
!!            open(unit=4,file='output/vtar100.dat')          
!!            open(unit=8,file='output/press_tar100.dat')          
!!            open(unit=9,file='output/vort_tar100.dat')     
c          else
c            open(unit=3,file='output/utar.dat')
c            open(unit=4,file='output/vtar.dat')
c            open(unit=8,file='output/press_tar.dat')
c            open(unit=9,file='output/vort_tar.dat')
c          endif
c        else
c          if( ifmm .eq. 1) then   
c            open(unit=3,file='output/utar_fmm.dat')      
c            open(unit=4,file='output/vtar_fmm.dat')      
c            open(unit=8,file='output/press_tar_fmm.dat')      
c            open(unit=9,file='output/vort_tar_fmm.dat') 
c          else   
c            open(unit=3,file='output/trap_utar.dat')      
c            open(unit=4,file='output/trap_vtar.dat')      
c            open(unit=8,file='output/trap_press_tar.dat')      
c            open(unit=9,file='output/trap_vort_tar.dat')
c          endif
c        endif
        
        do k = 1,nx*ny
c          write(3,1000) utar(k)
c          write(4,1000) vtar(k)
          write(8,1000) press_tar(k)
c          write(9,1000) vort_tar(k)
        enddo
c        close(unit=3)
c        close(unit=4)
        close(unit=8)
c        close(unit=9)
      enddo

 1000 format(E25.16)
 1010 format(ES20.4)
 1020 format(I4)



      end

