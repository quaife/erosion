      subroutine semiLagrangian(dt,ninner,nbodies,nouter,
     $    ifmm,ibary,x,y,den,ntrgets,xtar,ytar,xdep,ydep)
      implicit real*8 (a-h,o-z)

c     INPUT VARIABLES
      dimension x(ninner*nbodies),y(ninner*nbodies)
      dimension den(2*nouter + 2*ninner*nbodies + 3*nbodies)
      dimension xtar(ntargets),ytar(ntargets)

c     OUTPUT VARIABLES
      dimension xdep(ntargets),ydep(ntargets)

c     LOCAL VARIABLES
      press_tar(ntargets),vort_tar(ntargets)

c     Compute veocities
      call computerQoiTargets(ninner,nbodies,nouter,ibary,
     $    x,y,den,ntargets,xtar,ytar,utar,vtar,press_tar,vort_tar)

c     Compute first stage of midpoint rule
      do k = 1,ntargets
        xdep(k) = xtar(k) - 5.d-1*dt*utar
        ydep(k) = ytar(k) - 5.d-1*dt*vtar
      enddo

c     Compute veocities at the first stage
      call computerQoiTargets(ninner,nbodies,nouter,ibary,
     $    xdep,ydep,den,ntargets,xtar,ytar,utar,vtar,
     $    press_tar,vort_tar)
      do k = 1,ntargets
        xdep(k) = xdep(k) - dt*utar
        ydep(k) = ydep(k) - dt*vtar
      enddo


      return
      end
