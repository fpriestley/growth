program charge
  use constants_mod

  implicit none

  double precision :: dens,temp,efrac
  integer :: nsizes,zion
  double precision :: amin,amax,da
  double precision,allocatable :: agrain(:),qgrain(:)
  integer :: ncharge,zmax,zmin
  integer,allocatable :: zcharge(:)
  double precision,allocatable :: fcharge(:),j_pe(:),j_h(:),j_elec(:)
  character(len=100) :: nkfile
  integer :: nwav
  double precision,allocatable :: lambda(:),hv(:),nrad(:),krad(:),qabs(:),jrad(:)
  double precision,allocatable :: yield(:),detach(:),integrand(:)
  double precision :: wf,ebg
  integer :: getzmin,getzmax
  double precision :: photoyield,photdetach,jmat
  double precision :: integrate,limpe,limpd,lampet,lampdt
  double precision :: stick,jbar,coulombfac
  integer :: i,j,k
  

  dens = 30.
  temp = 100.
  efrac = 1.5e-4
  zion = 1

  nsizes = 100
  amin = 1e-3
  amax = 0.5

  wf = 8*ev
  ebg = 5*ev

  allocate(agrain(nsizes),qgrain(nsizes))

  da = (amax/amin)**(1./real(nsizes-1))

  do i=1,nsizes
     agrain(i) = amin * da**(i-1)
  end do

  nkfile = 'full_pyrmg100.lnk'

  call countlines(nkfile,nwav)

  allocate(lambda(nwav),hv(nwav),nrad(nwav),krad(nwav),qabs(nwav),jrad(nwav))
  allocate(yield(nwav),detach(nwav),integrand(nwav))

  call readnk(nkfile,nwav,lambda,nrad,krad)

  hv = hc/(1e-4*lambda)

  do i=1,nwav
     jrad(i) = jmat(lambda(i))
  end do

  do i=1,nsizes
     ! get Z range, setup arrays
     zmax = getzmax(agrain(i),wf)
     zmin = getzmin(agrain(i))
     ncharge = zmax-zmin+1
     allocate(zcharge(ncharge),fcharge(ncharge),j_pe(ncharge),j_h(ncharge),j_elec(ncharge))
     do j=1,ncharge
        zcharge(j) = zmin+j-1
     end do
     ! calculate Qabs
     call calcmie(nwav,lambda,nrad,krad,agrain(i),qabs)
     do j=1,ncharge
        ! calculate yield/photodetachment xsec
        yield = 0.
        detach = 0.
        limpe = lampet(agrain(i),zcharge(j),wf)
        limpd = lampdt(agrain(i),zcharge(j),wf,ebg)
        do k=1,nwav
           if (lambda(k) .le. limpe) yield(k) = photoyield(hv(k),zcharge(j),agrain(i),wf,ebg,krad(k))
           if ((zcharge(j) .lt. 0) .and. (lambda(k) .le. limpd)) detach(k) = photdetach(hv(k),zcharge(j),agrain(i),wf,ebg)
        end do
        ! calculate photoelectric rate
        integrand = fourpi*(1e-4*agrain(i))**2 * yield*qabs*1e-4*jrad/hv
        j_pe(j) = integrate(integrand,lambda,nwav,lambda(1),limpe)
        ! calculate photodetachment rate if -ve
        if (zcharge(j) .lt. 0) then
           integrand = detach*4e-4*jrad/hv
           j_pe(j) = j_pe(j) + integrate(integrand,lambda,nwav,lambda(1),limpd)
        end if
        ! calculate electron/hydrogen ion attachment
        j_elec(j) = efrac*dens*stick(agrain(i),zcharge(j),zmin) * sqrt(8*k_b*temp/(pi*m_e)) * pi*(1e-4*agrain(i))**2 &
             & * jbar(agrain(i),zcharge(j),-1,temp)
        j_h(j) = efrac*dens*sqrt(8*k_b*temp/(pi*m_p)) * pi*(1e-4*agrain(i))**2 * jbar(agrain(i),zcharge(j),1,temp)
     end do
     ! solve equilibrium charge distribution
     fcharge = 0.
     do j=1,ncharge
        if (ncharge .eq. 1) then
           fcharge(1) = 1.
           exit
        end if
        if (j_elec(j+1) .gt. 0.) then
           fcharge(j) = 1.
           do k=j+1,ncharge
              fcharge(k) = fcharge(k-1) * (j_h(k-1) + j_pe(k-1))/j_elec(k)
           end do
           exit
        end if
     end do
     fcharge = fcharge / sum(fcharge)
     ! calculate Coulomb focusing factor
     qgrain(i) = coulombfac(agrain(i),ncharge,zcharge,fcharge,temp,zion)
     write(*,'(ES10.3,2X,ES10.3)') agrain(i),qgrain(i)
     deallocate(zcharge,fcharge,j_pe,j_h,j_elec)
  end do

  deallocate(agrain,qgrain)
  deallocate(lambda,hv,nrad,krad,qabs,jrad)
  deallocate(yield,detach,integrand)

end program charge
