!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2017)                       !
!                                                     !
!=====================================================!
!                                                     !
!              SPHYNX: init_scenario.f90              !
!                                                     !
! Initial values for several model variables.         !
!=====================================================!

SUBROUTINE init_scenario

  USE parameters
  USE eosmodule

  IMPLICIT NONE
  INTEGER i,j,ii,jj
  DOUBLE PRECISION dmax,energytot,width,ener0,Aorb,v1,v2,m1,m2,vin
  DOUBLE PRECISION v3,v4
  INTEGER keytemp,keyerr
  DOUBLE PRECISION xrho,xye,xtemp,xtemp2
  DOUBLE PRECISION xenr,xprs,xent,xcs2,xdedt,xmunu
  DOUBLE PRECISION xdpderho,xdpdrhoe

  write(*,*) 'start initial scenario'
!--------------------  User can change this  --------------------

  
!----------------------------------------------------------------

  call calculate_hpowers
  !Desplacement and inicialization
  dmax=maxval(radius)                 !Estimation system size.
  rad=2.d0*sqrt(7.d0)*dmax            !Size of the tree root.
  despl=rad

  print *,'Rad:',rad

  a(:)=a(:)+despl

  call indexx(n,radius,indx)       !Creates sorted vector for first iteration

  !Initial velocity is zero?
  if(inivel0)v=0.d0
  
  
#ifdef EQMASS
     TD=1.d0/sqrt(g*dble(n)*masspart/dmax)
#else
     TD=1.d0/sqrt(g*sum(mass)/dmax)
#endif
 

  print *,'Rad:',rad

  call indexx(n,radius,indx)       !Creates sorted vector for first iteration



  !Volume Elements

#ifdef EQMASS
  !Set to the STD until 5 iteration before NR update starts
  if(iterini.le.liniNR-5) then
     xmass(:)=masspart
  else
     xmass(:)=masspart/promro(:)
  endif
#else
  !Set to the STD until 5 iteration before NR update starts
  if(iterini.le.liniNR-5) then
     xmass(:)=mass(:)
  else
     xmass(:)=mass(:)/promro(:)
  endif
#endif

  write(*,*) 'End initial scenario'

  RETURN
END SUBROUTINE init_scenario
