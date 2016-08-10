      program stokesDriver
      implicit real*8 (a-h,o-z)

      parameter (ninner = 256)

      dimension x(ninner),y(ninner)
      dimension shear_stress(ninner)


      twopi = 8.d0*datan(1.d0)
      dtheta = twopi/dble(ninner)
      do k = 1,ninner
        theta = dble(k-1)*dtheta
        x(k) = 5.d-1*dcos(theta)
        y(k) = 5.d-1*dsin(theta)
c        rad = 2.d-1 + 8.d-2*dcos(5*theta)
c        x(k) = rad*dcos(theta)
c        y(k) = rad*dsin(theta)
      enddo
c     parameterize from left most point and proceed up the top of the
c     curve, throught the back point, and back to the leading
c     singularity point

      call stokesSolver(ninner,x,y,shear_stress)
c     pass in number of points and x and y coordinates and return the
c     shear stress on the boundary

      end

