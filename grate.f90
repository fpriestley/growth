double precision function growthrate(n_gas,t_gas,rho,m_elem,f_elem,gfrac,fgas)
  use constants_mod

  implicit none

  double precision,intent(in) :: n_gas,t_gas,rho,m_elem,f_elem,gfrac,fgas
  double precision :: vx,n_elem,stick

  vx = sqrt(k_b*t_gas/m_elem/twopi)
  n_elem = n_gas * f_elem * fgas

  growthrate = 1e4*stick(t_gas)*n_elem*vx*m_elem/rho/gfrac/4.

end function growthrate

double precision function stick(t)
  use constants_mod

  implicit none

  double precision,intent(in) :: t
  double precision :: tmax

  tmax = 300.

  if (t .gt. tmax) then
     stick = 0.
  else
     stick = 1.
  end if

end function stick

subroutine growsizes(nsizes,ngrain,qgrain,agrain,dt,n_gas,t_gas,rho,m_elem,f_elem,gfrac,fgas)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: qgrain(nsizes),agrain(nsizes),dt,n_gas,t_gas,rho,m_elem,f_elem,gfrac,fgas
  double precision,intent(inout) :: ngrain(nsizes)
  double precision :: newngrain(nsizes),anew(nsizes)
  integer :: i,j
  double precision :: growthrate,dadt,gr,a1,a2,da,fa1,fa2

  gr = growthrate(n_gas,t_gas,rho,m_elem,f_elem,gfrac,fgas)

  newngrain = 0.
  anew = 0.

  do i=1,nsizes-1
     if (charge) then
        dadt = qgrain(i) * gr
     else
        dadt = gr
     end if
     anew(i) = agrain(i) + dt*dadt
     do j=i+1,nsizes
        if (agrain(j) .ge. anew(i)) then
           a1 = agrain(j-1)
           a2 = agrain(j)
           da = a2-a1
           fa1 = 1.-(anew(i)-a1)/da
           fa2 = 1.-(a2-anew(i))/da
           newngrain(j) = newngrain(j) + fa2*ngrain(i)
           newngrain(j-1) = newngrain(j-1) + fa1*ngrain(i)
           exit
        end if
     end do
  end do

  newngrain(nsizes) = newngrain(nsizes) + ngrain(nsizes)

  ngrain = newngrain

end subroutine growsizes

subroutine destroy(nsizes,ngrain,agrain,dt,tdest,anorm)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),dt,tdest,anorm
  double precision,intent(inout) :: ngrain(nsizes)
  double precision :: ta
  integer :: i

  do i=1,nsizes
     ta = tdest !* (agrain(i)/anorm)
     ngrain(i) = ngrain(i)*(1 - dt/ta)
  end do

end subroutine destroy

double precision function calcgtimescale(n_gas,t_gas,rho,m_elem,f_elem,gfrac,avg_a3)
  use constants_mod

  implicit none

  double precision,intent(in) :: n_gas,t_gas,rho,m_elem,f_elem,gfrac,avg_a3
  double precision :: vx,stick

  vx = sqrt(k_b*t_gas/m_elem/twopi)

  calcgtimescale = 3. * stick(t_gas) * f_elem * m_elem * vx * n_gas
  calcgtimescale = calcgtimescale / (rho * gfrac * 1e-4*avg_a3)
  calcgtimescale = 1./calcgtimescale

end function calcgtimescale
