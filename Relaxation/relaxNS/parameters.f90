!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2018)                       !
!                                                     !
!=====================================================!
!                                                     !
!                SPHYNX: parameters.f90               !
!                                                     !
! Main module with all parameters and global variables!
!                                                     !
!=====================================================!

MODULE parameters

  IMPLICIT NONE
  DOUBLE PRECISION rad,ddmmax,rocentral,despl
  DOUBLE PRECISION dt,dtnma,tt,dtold
  DOUBLE PRECISION mtotal,cmax
  DOUBLE PRECISION ad1,ad2,ad3,adw1,adw2,adw3,pk4
  DOUBLE PRECISION refrad,refax,refay,refaz
  DOUBLE PRECISION maxtstepaux,timeinit,timefin,masspart
  DOUBLE PRECISION masscloudinic,TD

  INTEGER l,n,npini,npend,id,iti,nproc
  INTEGER ntotal,lastcoldim,new_comm

  LOGICAL escribe,findballmass
  LOGICAL checkdens,initLSvar

  INTEGER,PARAMETER::n1 = 100000    !Number of particles for star 1 (equal for star 2)
  INTEGER,PARAMETER::nmax = n1    !Maximum number of particles

  INTEGER,PARAMETER::n2=2*nmax       !Multiples of nmax
  INTEGER,PARAMETER::n3=3*nmax
  INTEGER,PARAMETER::n4=4*nmax
  INTEGER,PARAMETER::n5=5*nmax
  INTEGER,PARAMETER::n7=7*nmax

  INTEGER,PARAMETER::nodes=7         !Amount of nodes.
  INTEGER,PARAMETER::dim = 3         !Spatial imensions
  INTEGER,PARAMETER::nnl = 15000     !Iterations
  INTEGER,PARAMETER::iodelay = 100   !Steps until output
  INTEGER,PARAMETER::nt = 3*nmax     !Max # of nodes in tree
  INTEGER,PARAMETER::ntree = nmax/2  !Max depth in tree
  INTEGER,PARAMETER::ncubes = 1      !Cubes for periodic BC
                                     !1 - No periodic BC
                                     !9 - Full periodic BC 2D
                                     !27- Full periodic BC 3D
  INTEGER,PARAMETER::kernel = 3      !Choosing kernel
                                     !1 - Cubic Spline
                                     !2 - Harmonic
                                     !3 - Harmonic interpolated
  LOGICAL,PARAMETER::oldnorm=.false. ! Choose old interpolation for norm
  INTEGER,PARAMETER::NRitermax = 10  !Maximum NR iterations for h
  INTEGER,PARAMETER::liniNR = nnl+10     !Iteration to start NR for h

  DOUBLE PRECISION,PARAMETER::scale=1.d0      !spatial scale (cm)
  DOUBLE PRECISION,PARAMETER::xbox=.25d0      !Box size x
  DOUBLE PRECISION,PARAMETER::ybox=.25d0      !Box size y
  DOUBLE PRECISION,PARAMETER::zbox=1.d0       !Box size z
  DOUBLE PRECISION,PARAMETER::ni0=100.d0      !Number of neighbors
  DOUBLE PRECISION,PARAMETER::nvmax=160.d0    !Maximum number of neighbors(Int.)
  DOUBLE PRECISION,PARAMETER::nvmin=60.d0     !Minimum number of neighbors(Int.)
  DOUBLE PRECISION,PARAMETER::fixfactor=3.d0  !Increased % for the neigh. array
  DOUBLE PRECISION,PARAMETER::hfac=1023.d0      !Factor for h update
  DOUBLE PRECISION,PARAMETER::hexp=1.d0/10.d0  !Exponent for h update
  DOUBLE PRECISION,PARAMETER::tol = 0.5d0     !Tolerance for opening nodes
  DOUBLE PRECISION,PARAMETER::dkcour = 0.2d0  !Factor for courant timestep
  DOUBLE PRECISION,PARAMETER::dkro = .06d0    !Factor for density timestep
  DOUBLE PRECISION,PARAMETER::dkacc = .02d0   !Factor for acceleration timestep
  DOUBLE PRECISION,PARAMETER::deltahindex=0.d0!Maximum index change
  DOUBLE PRECISION,PARAMETER::lambdac=0.5d0   !Parameter for index change
  DOUBLE PRECISION,PARAMETER::nbaseline=6.d0  !Minimum index
  DOUBLE PRECISION,PARAMETER::alfamax=2.d0    !Artificial viscosity constant
  DOUBLE PRECISION,PARAMETER::alfau=0.d0      !Artificial viscosity thermal. Put it to zero if not wanted
  DOUBLE PRECISION,PARAMETER::alfamin=0.05d0  !AV floor (switches)
  DOUBLE PRECISION,PARAMETER::betaAV=3.d0     !Artificial viscosity constant beta
  DOUBLE PRECISION,PARAMETER::decay_constant=0.2     !Alfa decay (switches)
  DOUBLE PRECISION,PARAMETER::p0=1.d0         !Initial constant pressure
  !Min. Atwood number in ramp function in momentum equation (crossed/uncrossed selection)
  !Complete uncrossed option  (Atmin>=1.d50, Atmax it doesn't matter).
  !Complete crossed (Atmin and Atmax negative)
  DOUBLE PRECISION,PARAMETER::Atmin=0.1d0
  DOUBLE PRECISION,PARAMETER::Atmax=0.2d0
  DOUBLE PRECISION,PARAMETER::ramp=1.d0/(Atmax-Atmin)

  DOUBLE PRECISION,PARAMETER::maxtstep=1.d-3  !Biggest allowed timestep

  LOGICAL,PARAMETER::flags = .true.       !Activate verbosity
  LOGICAL,PARAMETER::relax = .false.      !Activate angular relaxation
  LOGICAL,PARAMETER::gravity = .true.     !Calculate gravity
  LOGICAL,PARAMETER::estabil = .true.     !Write central density and radius
  LOGICAL,PARAMETER::conserv = .true.     !Tracking of conservation laws
  LOGICAL,PARAMETER::timesave = .true.    !Tracking of timesteps
  LOGICAL,PARAMETER::balsara = .false.    !Balsara corrections
  LOGICAL,PARAMETER::switches = .false.   !AV switches
  LOGICAL,PARAMETER::gradh = .true.       !Includes grad-h terms
  LOGICAL,PARAMETER::conduction=.false.   !Thermal conduction (not in connection with AV)
  LOGICAL,PARAMETER::std_VE=.true.        !Use xmass=mass(i) or new  xmass=mass(i)/ro0(i)
  LOGICAL,PARAMETER::outini=.true.        !Write output file at first iteration
  LOGICAL,PARAMETER::timecour = .true.    !Courant timestep
  LOGICAL,PARAMETER::timero = .false.     !Densidad timestep
  LOGICAL,PARAMETER::timeacc = .false.    !Acceleration timestep
  LOGICAL,PARAMETER::KH=.false.           !set KH test in EOS
  LOGICAL,PARAMETER::sortneigh=.false.    !Sort list of neighbors
  LOGICAL,PARAMETER::gwrad=.false.        !Calculate gravitational wave emission
  LOGICAL,PARAMETER::initLS = .true.      !Initialize U using tinterp

  CHARACTER*43,PARAMETER::formatin='(37(1x,es23.16))' !Input format
  CHARACTER*43,PARAMETER::formatout='(1x,i7,45(1x,es12.5),1x,i2,1x,l1)'  !Output format

  !RESTART parameters
  INTEGER,PARAMETER::iterini = 1              !Begining of first iteration
  INTEGER,PARAMETER::delay = 0                !timesteps of delay to h adjustment
  DOUBLE PRECISION,PARAMETER::timeini = 0.d0  ! Origin of time
  DOUBLE PRECISION,PARAMETER::deltaini = 1.d-9    !Initial timestep
  LOGICAL,PARAMETER::inivel0 = .true.             !Initial velocity=0
  LOGICAL,PARAMETER::inienergy = .false.          !Initialize U with EOS
  CHARACTER*20,PARAMETER::inputfile='neutron_star'  !Input file
  CHARACTER*100,PARAMETER::comments=&
       &'Neutron star relaxation'
  !EOS
  DOUBLE PRECISION, PARAMETER::gamma=5.d0/3.d0
  INTEGER,PARAMETER::eos = 4 !Choose EOS
                              !1 - Ideal gas with gamma with Temp
                              !2 - Ideal gas with gamma with Eint
                              !3 - Polytropic
                              !4 - tabulated EoS (Driver by O'Connor et al.)
                              !5 - piecewise polytrope
                              !6 - LS220 (Driver by Matthias Liebendörfer)

  CHARACTER*24,PARAMETER::EOStab='../eos/LS220.h5'      !Tabulated EOS
  DOUBLE PRECISION,PARAMETER::kpol = 4.3130893169156435d-2 !Polytropic constant
  DOUBLE PRECISION,PARAMETER::gammapol = 1.d0+1.d0/0.7d0           !Polytropic exponent
  
  !piecewise polytrope parameters
  DOUBLE PRECISION,PARAMETER::kpol0= 3.59389d13
  DOUBLE PRECISION,PARAMETER::gammapol0=1.3569d0
  DOUBLE PRECISION,PARAMETER::gammapol1=3.446d0
  DOUBLE PRECISION,PARAMETER::gammapol2=3.572d0
  DOUBLE PRECISION,PARAMETER::gammapol3=2.887d0
  DOUBLE PRECISION,PARAMETER::press1=10**34.495
  DOUBLE PRECISION,PARAMETER::rhoth1=10.d0**14.7d0
  DOUBLE PRECISION,PARAMETER::rhoth2=1.000000d15
  DOUBLE PRECISION,PARAMETER::kpol1= press1/(rhoth1**gammapol1)
  DOUBLE PRECISION,PARAMETER::rhoth0=(kpol1/kpol0)**(1.d0/(gammapol0-gammapol1))
  DOUBLE PRECISION,PARAMETER::kpol2= kpol1*rhoth1**(gammapol1-gammapol2)
  DOUBLE PRECISION,PARAMETER::kpol3= kpol2*rhoth2**(gammapol2-gammapol3)
  DOUBLE PRECISION,PARAMETER::pth0=kpol1*rhoth0**gammapol1
  DOUBLE PRECISION,PARAMETER::pth1=kpol2*rhoth1**gammapol2
  DOUBLE PRECISION,PARAMETER::pth2=kpol3*rhoth2**gammapol3

  !Physical and math constants
  DOUBLE PRECISION,PARAMETER::g=6.6726d-8     !Gravitation constant
  DOUBLE PRECISION,PARAMETER::cvel=2.997924562d10       !Speed of light
  DOUBLE PRECISION,PARAMETER::pi=3.14159265358979323846d0
  DOUBLE PRECISION,PARAMETER::pihalf=pi/2.d0
  DOUBLE PRECISION,PARAMETER::pim1=1.d0/pi
  DOUBLE PRECISION,PARAMETER::third=1.d0/3.d0
  DOUBLE PRECISION,PARAMETER::twothird=2.d0/3.d0
  DOUBLE PRECISION,PARAMETER::msun=1.989d33
  DOUBLE PRECISION,PARAMETER::rsun=6.9598d10
  DOUBLE PRECISION,PARAMETER::MeV=1.602192d-6
  DOUBLE PRECISION,PARAMETER::mb=1.674d-24    !Baryon mass
  DOUBLE PRECISION,PARAMETER::K2MeV=8.617342857d-11 !Conversion from K to MeV
  DOUBLE PRECISION,PARAMETER::rgasid=8.317d7
  DOUBLE PRECISION,PARAMETER::MASSN=939.5731d0!Neutron mass in MeV for LSEOS
  DOUBLE PRECISION,PARAMETER::ZE=8.9d0        !Zero energy for LSEOS

  !Arrays
  INTEGER,DIMENSION(nmax)::indx,nvi
  INTEGER,DIMENSION(0:nt)::nh1,nh2
  INTEGER,DIMENSION(nt)::nnod,npa
  INTEGER,ALLOCATABLE,DIMENSION(:,:)::neighbors,cube

#ifdef EQMASS
#else
  DOUBLE PRECISION,DIMENSION(nmax)::mass
#endif
  DOUBLE PRECISION,DIMENSION(nmax)::temp,promro,ro2,u,p,c,pro,promro0
  DOUBLE PRECISION,DIMENSION(nmax)::mue,mui,cv,dudv,dpdt,ro0,sumwhro0
  DOUBLE PRECISION,DIMENSION(nmax)::yelec,s,alfa
  DOUBLE PRECISION,DIMENSION(nmax)::mindist,avgro,xmass,sumkx
  DOUBLE PRECISION,DIMENSION(nmax)::h,hd,h2,h3,hm1,hm2,hm3,hm4,hm5,omega
  DOUBLE PRECISION,DIMENSION(nmax)::energy,avisc,radius,sumwh,ballmass
  DOUBLE PRECISION,DIMENSION(nmax)::vol
  DOUBLE PRECISION,DIMENSION(nmax)::energyold,aviscold,aviscu,aviscuold
  DOUBLE PRECISION,DIMENSION(nmax)::acmod
  DOUBLE PRECISION,DIMENSION(nmax)::c11,c12,c13,c22,c23,c33
  DOUBLE PRECISION,DIMENSION(nmax)::magcorr,gradcorr,maxvsignal
  DOUBLE PRECISION,DIMENSION(nmax)::kappa,kappa0,vturb
  DOUBLE PRECISION,DIMENSION(nmax)::indice,pk,dpk,mark_ramp
  DOUBLE PRECISION,DIMENSION(nmax)::checkInorm,checkdeltax,checkdeltay,checkdeltaz

  DOUBLE PRECISION,DIMENSION(20001)::uk,fk,fdk1,fdk2
  DOUBLE PRECISION,dimension(25)::timepass

  DOUBLE PRECISION,DIMENSION(nmax)::diff

  DOUBLE PRECISION,DIMENSION(nt)::xce,yce,zce,dxx,dxx2,xmc,ymc,zmc
  DOUBLE PRECISION,DIMENSION(nt)::qxx,qxy,qxz,qyy,qyz,qzz,trq,mcm

  DOUBLE PRECISION,DIMENSION(nmax*dim)::a,a0,v,v0,f,ac,gradp

  DOUBLE PRECISION,DIMENSION(dim)::cm

  DOUBLE PRECISION,DIMENSION(nmax)::ugrav,divv,curlv

  LOGICAL,DIMENSION(nmax)::wind,ready,mark
  LOGICAL,DIMENSION(nodes)::localready
  
  DOUBLE PRECISION::sepa,dEdt,initAorb
  DOUBLE PRECISION,DIMENSION(dim)::cm1,cm2

END MODULE parameters
