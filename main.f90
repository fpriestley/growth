program growth
  use constants_mod

  implicit none

  character(len=100) :: name,coulombname
  integer :: nsizes
  double precision :: amin,amax,mrn
  double precision,allocatable :: agrain(:),ngrain(:),awidth(:),qgrain(:)
  double precision :: depletion,rho_g,gfrac
  double precision :: m_elem,f_elem,fgas,fdust,mdust,fgas_ng
  double precision :: avg_a3,gtimescale,dfdt_ng,gtimescale0
  double precision :: calcavg_a3,calcgtimescale
  double precision :: n_gas,t_gas
  double precision :: time,time_max,dt
  double precision :: massingrains
  double precision :: vturb,shateff
  double precision :: tdest,anorm
  integer :: nq,i
  double precision,allocatable :: aq(:),fq(:)

  open(unit=10,file='input/input.dat',status='old')

  read(10,*) name
  read(10,*) n_gas
  read(10,*) t_gas
  read(10,*) tdest
  read(10,*) nsizes
  read(10,*) amin
  read(10,*) amax
  read(10,*) coulombname
  read(10,*) vturb
  read(10,*) shateff

  close(unit=10)

  call openfiles(name)

  mrn = -3.5
  rho_g = 3.13
  gfrac = 0.165
  f_elem = 3.24e-5
  m_elem = 28*m_p

  fgas = 10**(-0.5)
  fgas_ng = fgas
  fdust = 1. - fgas

  tdest = tdest * 1e6*yr
  anorm = 0.01

  time_max = 5e8*yr
  dt = 1e4*yr
  time = 0.

  allocate(agrain(nsizes))
  allocate(ngrain(nsizes))
  allocate(awidth(nsizes))
  allocate(qgrain(nsizes))

  call countlines('input/'//trim(coulombname),nq)

  allocate(aq(nq))
  allocate(fq(nq))

  call readfactor(coulombname,nq,aq,fq)
  call makesizedist(nsizes,agrain,ngrain,awidth,qgrain,amin,amax,mrn,nq,aq,fq)
  call initdepletion(nsizes,agrain,ngrain,rho_g,f_elem,m_elem,gfrac,fdust,n_gas)
  
  deallocate(aq)
  deallocate(fq)

  avg_a3 = calcavg_a3(nsizes,agrain,ngrain,qgrain)
  gtimescale = calcgtimescale(n_gas,t_gas,rho_g,m_elem,f_elem,gfrac,avg_a3)
  gtimescale0 = gtimescale
  mdust = massingrains(nsizes,agrain,ngrain,rho_g)

  call writedepletion(time/yr,fgas,fgas_ng,mdust,gtimescale/yr)
  call writesizes(nsizes,ngrain)

  do
     time = time + dt
     if (time .gt. time_max) exit
     call growsizes(nsizes,ngrain,qgrain,agrain,dt,n_gas,t_gas,rho_g,m_elem,f_elem,gfrac,fgas)
     if (lgdestroy) call destroy(nsizes,ngrain,agrain,dt,tdest,anorm)
     if (lgshatter) call shatter(nsizes,agrain,ngrain,dt,vturb,rho_g,shateff)
     call updatedepletion(nsizes,agrain,ngrain,rho_g,f_elem,m_elem,gfrac,fgas,fdust,n_gas)
     mdust = massingrains(nsizes,agrain,ngrain,rho_g)
     fgas_ng = fgas_ng - dt*fgas_ng*(1.-fgas_ng)/gtimescale0
     if (lgdestroy) fgas_ng = fgas_ng + dt*(1.-fgas_ng)/tdest
     avg_a3 = calcavg_a3(nsizes,agrain,ngrain,qgrain)
     gtimescale = calcgtimescale(n_gas,t_gas,rho_g,m_elem,f_elem,gfrac,avg_a3)
     call writedepletion(time/yr,fgas,fgas_ng,mdust,gtimescale/yr)
     call writesizes(nsizes,ngrain)
     if (fgas .le. 0.0d0) stop
  end do

  deallocate(agrain)
  deallocate(ngrain)
  deallocate(awidth)
  deallocate(qgrain)

  call closefiles()

end program growth
