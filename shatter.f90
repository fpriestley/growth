subroutine shatter(nsizes,agrain,ngrain,dt,vturb,rho,shateff)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),dt,vturb,rho,shateff
  double precision,intent(inout) :: ngrain(nsizes)
  double precision :: newn(nsizes),tempn(nsizes),area,shatfrac,shatrate,shatmass
  double precision :: integrate,norm
  integer :: i,j

  newn = 0.

  do i=1,nsizes
     area = 0.
     do j=1,nsizes
        area = area + ngrain(j)*(1e-4*agrain(i) + 1e-4*agrain(j))**2
     end do
     shatrate = shateff*area*pi*vturb
     shatfrac = min(shatrate*dt,1.)
     shatmass = shatfrac*ngrain(i)*fourthirdpi*rho*(1e-4*agrain(i))**3
     newn(i) = newn(i) + (1-shatfrac)*ngrain(i)
     tempn = 0.
     do j=1,i
        tempn(j) = (1e-4*agrain(j))**(-2.5)
     end do
     norm = shatmass/sum(fourthirdpi*rho*tempn*(1e-4*agrain)**3)
     tempn = tempn*norm
     newn = newn + tempn
  end do

  ngrain = newn

end subroutine shatter
