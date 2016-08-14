module mymod 

contains
   function foo(xx)
      real*8 :: foo, xx
      foo = 2*xx
   end function foo

   subroutine add(xx,yy)
      real*8, intent(in) :: xx(2)
      real*8, intent(out) :: yy(2)
      yy(1) = xx(1)+xx(2)
      yy(2) = xx(1)*xx(2)
   end subroutine add

end module mymod