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
  !define parameters for seperation and initial orbital velocity
  dmax=maxval(radius)
  print *,'initial radius:',dmax
  Aorb=4.d0*dmax
  initAorb=Aorb
  m1=dble(n1)*masspart
  m2=dble(n-n1)*masspart
  v1=sqrt(g*(m1)/(2.d0*Aorb)) !orbital velocity
  !vin=v1/29979245800.d0
  !print *,'orbital velocity [c]:',vin
  v2=-v1
  v3=64.d0/5.d0*1/(cvel**5)*(g*m1/Aorb)**3.d0
  v4=-v3
  
  !seperate star 1 and star 2 by Aorb
  !open(10,file='stars.txt')
     do i=1,n1
       j=i+n1
       ii=1+dim*(i-1)
       jj=1+dim*(j-1)
       a(ii)=a(ii)-(Aorb/2.d0)
       a(jj)=a(jj)+(Aorb/2.d0)
       v(ii)=v3
       v(jj)=v4
       v(ii+1)=v1
       v(jj+1)=v2
  !     !write(10,*) a(ii),a(1+ii),a(2+ii),a(jj),a(jj+1),a(jj+2)
     enddo
  !close(10)
  
  !find new masscenter
  call masscenter(1,n)
  print *,'masscenter:', cm(:)
  do i=1,n
      ii=1+dim*(i-1)
      a(ii)=a(ii)-cm(1)
      a(ii+1)=a(ii+1)-cm(2)
      a(ii+2)=a(ii+2)-cm(3)
      radius(i)=sqrt(a(ii)*a(ii)+a(ii+1)*a(ii+1)+a(ii+2)*a(ii+2))
  enddo

  !Desplacement and inicialization
  dmax=maxval(radius)                 !Estimation system size.
  print *,'system size:',dmax
  rad=2.d0*sqrt(7.d0)*dmax            !Size of the tree root.
  despl=rad
  a(:)=a(:)+rad

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
