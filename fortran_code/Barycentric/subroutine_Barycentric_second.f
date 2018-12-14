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
     $    ntar,xtar,ytar,utar,dutar,ddutar)
     
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
     $    ntar,xtar,ytar,utar,dutar,ddutar)     

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
     $    ntar,xtar,ytar,utar,dutar,ddutar)     

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
     $    ntar,xtar,ytar,utar,dutar,ddutar)

      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(dutar(j))
        u1ytar(j) = u1ytar(j) - dimag(dutar(j))
        endif
      enddo

      call compute_bdval_in(npts,xsou,ysou,tau2,bdval)

      call laplaceInteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,dutar,ddutar)

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
     $      ntar,xtar,ytar,utar,dutar,ddutar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     $             ddutar(ntar)

      complex *16 eye,czero
      complex*16 zsou(npts),dzsou(npts)
      complex *16 num,den,num1
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

      do j = 1,ntar
        ztar = xtar(j) + eye*ytar(j)
        num = czero
        den = czero
        do k = 1,npts
          num = num + bdval(k)/(zsou(k) - ztar)*dzsou(k)
          den = den + 1.d0/(zsou(k) - ztar)*dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.5) of
c       of Barnett, Wu, Veerapaneni

        utar(j) = num/den


c       Start computing derivative (equation (3.10)).  Term in
c       denominator is unchanged
        num = czero
        do k = 1,npts
          num = num + (bdval(k) - utar(j))/
     $          (zsou(k) - ztar)**2.d0*dzsou(k)
        enddo
c       compute numerator and denominator term in equation (3.5) of
c       Barnett, Wu, Veerapaneni
        dutar(j) = num/den

c      Computing the second derivative of DLP for the stress tensor

       num = czero
       num1 = czero
       
         Do k=1,npts       
             num1 = num1+(bdval(k)-utar(j)-dutar(j)*(zsou(k)-ztar))/
     $               (zsou(k)-ztar)**3.d0*dzsou(k)
         enddo
         ddutar(j)=2.d0*num1/den
       enddo

      end

c****************************************************************
      subroutine StokesExteriorDLP(npts,nptsc,nbeta,xsou,ysou,
     &          densx,densy,px,py,ntar,xtar,ytar,isou,n,nder,
     &          utar1x,utar1y,u1xtar,u1ytar,u2xtar,u2ytar)!,press_tar)
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
        tau1(j) = (densx(j) + eye*densy(j))*px(j)/outnor
        tau2(j) = (densx(j) + eye*densy(j))*py(j)/outnor
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
     $    ntar,xtar,ytar,utar,dutar,ddutar)
     
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
     $    ntar,xtar,ytar,utar,dutar,ddutar)     

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
     $    ntar,xtar,ytar,utar,dutar,ddutar)     

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
     $    ntar,xtar,ytar,utar,dutar,ddutar)
      
      do j = 1,ntar
        utar1x(j) = utar1x(j) + dreal(utar(j))
        
        if(nder .eq. 2) then        
        u1xtar(j) = u1xtar(j) + dreal(dutar(j))
        u1ytar(j) = u1ytar(j) - dimag(dutar(j))
        endif
        
      enddo

      call compute_bdval_ex(npts,xsou,ysou,tau2,bdval)

      call laplaceExteriorHolomorphic(npts,xsou,ysou,bdval,
     $    ntar,xtar,ytar,utar,dutar,ddutar)
      
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
     $      ntar,xtar,ytar,utar,dutar,ddutar)
      implicit real*8 (a-h,o-z)

      dimension xsou(npts),ysou(npts)
      dimension xtar(ntar),ytar(ntar)
      complex *16 bdval(npts)
      complex *16 utar(ntar),dutar(ntar),
     &            ddutar(ntar) 

      complex *16 eye,czero,a
      complex*16 zsou(npts),dzsou(npts)
      complex *16 num,den,num1
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
          num = num + bdval(k)/(zsou(k) - ztar)*dzsou(k)
c          den = den + 1.d0/((zsou(k) - dcmplx(a))
c     $          *(zsou(k) - ztar))*dzsou(k)
          den = den + 1.d0/(zsou(k) - ztar)*dzsou(k)
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
        if( CDABS(zsou(k)-ztar) .gt. 1.d-10) then
          num = num + (bdval(k) - utar(j))/
     $          (zsou(k) - ztar)**2.d0*dzsou(k)
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
        
c      Computing the second derivative of DLP for the stress tensor
        
       num = czero
       num1 = czero
       
         Do k=1,npts       
             num1 = num1+(bdval(k)-utar(j)-dutar(j)*(zsou(k)-ztar))/
     $               (zsou(k)-ztar)**3.d0*dzsou(k)
         enddo
         
         ddutar(j)=2.d0*num1/(ztar - dcmplx(a))/den
       enddo

      end
     
