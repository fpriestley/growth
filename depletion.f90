double precision function massingrains(nsizes,agrain,ngrain,rho)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),ngrain(nsizes),rho
  double precision :: mass_a(nsizes)
  integer :: i

  do i=1,nsizes
     mass_a(i) = ngrain(i) * fourthirdpi * (1e-4*agrain(i))**3 * rho
  end do

  massingrains = sum(mass_a)

end function massingrains

subroutine initdepletion(nsizes,agrain,ngrain,rho,f_elem,m_elem,gfrac,fdust,ngas)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),rho,f_elem,m_elem,gfrac,fdust,ngas
  double precision,intent(inout) :: ngrain(nsizes)
  double precision :: massingrains,avbl_mass,startmass

  avbl_mass = fdust*f_elem*m_elem/gfrac*ngas
  startmass = massingrains(nsizes,agrain,ngrain,rho)

  ngrain = avbl_mass * ngrain/startmass

end subroutine initdepletion

subroutine updatedepletion(nsizes,agrain,ngrain,rho,f_elem,m_elem,gfrac,fgas,fdust,n_gas)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: agrain(nsizes),ngrain(nsizes),rho,f_elem,m_elem,gfrac,n_gas
  double precision,intent(inout) :: fdust,fgas
  double precision :: massingrains,dustmass,totmass

  dustmass = gfrac*massingrains(nsizes,agrain,ngrain,rho)
  totmass = f_elem*m_elem*n_gas

  fdust = dustmass/totmass
  fgas = (totmass-dustmass)/totmass

end subroutine updatedepletion
