subroutine interpolate(x,y,xtab,ytab,n)

  implicit none

  integer,intent(in) :: n
  double precision,intent(in) :: x,xtab(n),ytab(n)
  double precision,intent(out) :: y
  integer :: i
  double precision :: dx,grad

  if (x .lt. xtab(1)) then
     y = ytab(1)
  else if (x .gt. xtab(n)) then
     y = ytab(n)
  else
     do i=1,n-1
        if (xtab(i+1) .ge. x) then
           grad = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
           dx = x-xtab(i)
           y = ytab(i) + grad*dx
           exit
        end if
     end do
  end if

end subroutine interpolate
