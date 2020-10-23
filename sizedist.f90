subroutine makesizedist(nsizes,agrain,ngrain,awidth,qgrain,amin,amax,mrn,nq,aq,fq)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes,nq
  double precision,intent(in) :: amin,amax,mrn,aq(nq),fq(nq)
  double precision,intent(out) :: agrain(nsizes),ngrain(nsizes),awidth(nsizes),qgrain(nsizes)
  double precision :: atemp(0:nsizes),f(0:nsizes),da
  double precision :: integrate,norm
  integer :: i

  da = (amax/amin)**(1./real(nsizes))

  do i=0,nsizes
     atemp(i) = amin * da**i
     f(i) = atemp(i)**mrn
  end do

  norm = integrate(f,atemp,nsizes+1,atemp(0),atemp(nsizes))
  f = f/norm

  do i=1,nsizes
     agrain(i) = 0.5*(atemp(i)+atemp(i-1))
     ngrain(i) = 0.5*(f(i)+f(i-1))*(atemp(i)-atemp(i-1))
     awidth(i) = atemp(i)-atemp(i-1)
  end do

  call calcfactor(nsizes,agrain,qgrain,nq,aq,fq)

  write(32,"(1000ES10.3)") agrain

end subroutine makesizedist

double precision function calcavg_a3(nsizes,agrain,ngrain,qgrain)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),ngrain(nsizes),qgrain(nsizes)
  double precision :: a3,a2

  a3 = sum(ngrain*agrain**3)
  if (charge) then
     a2 = sum(ngrain*qgrain*agrain**2)
  else
     a2 = sum(ngrain*agrain**2)
  end if

  calcavg_a3 = a3/a2

end function calcavg_a3
