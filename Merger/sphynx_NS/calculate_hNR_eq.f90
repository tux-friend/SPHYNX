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
!             SPHYNX: calculate_hNR_eq.f90            !
!                                                     !
! Newton-Raphson for h and n.                         !
!=====================================================!

  SUBROUTINE calculate_hNR_eq
    
    USE parameters

    IMPLICIT NONE

    INTEGER i
    DOUBLE PRECISION lambda0,exponential1,exponential2,xplus,xminus,lambda1
    DOUBLE PRECISION g1,g11,g12,g2,g21,g22,det,detg1,detg2,deltah,deltan

    print *,'calculate h Newton-Rapshon',id
    pk=0.d0
    dpk=0.d0    
    dnro=0.d0

!$omp parallel private(i,lambda0,exponential1,exponential2,xplus,&
!$omp                  xminus,lambda1,g1,g11,g12,g2,g21,g22,det,&
!$omp                  detg1,detg2,deltah,deltan)
!$omp do schedule(static)
    do i=npini,npend
       roprom(i)=exp(roprom(i))
       lambda0=roprom(i)/promro(i)
       if(lambda0.ge.1.d0) then
          exponential1=exp((lambda0-1.d0)/lambdac)
          exponential2=1.d0/exponential1
          xplus=exponential1+exponential2
          xminus=exponential1-exponential2
          dnro(i)=deltahindex*2.d0*xminus/xplus**2/lambdac
          dnro(i)=dnro(i)*lambda0/promro(i)*(1.d0/(dble(nvi(i))+1.d0)-1.d0)
       else
          lambda1=1.d0/lambda0
          exponential1=exp((lambda1-1.d0)/lambdac)
          exponential2=1.d0/exponential1
          xplus=exponential1+exponential2
          xminus=exponential1-exponential2
          dnro(i)=deltahindex*2.d0*xminus/xplus**2/lambdac
          dnro(i)=dnro(i)*lambda1/promro(i)*(1.d0-1.d0/(dble(nvi(i))+1.d0))
       endif

       ! Find the zero of functions  g1 and g2 
       g1=ballmass(i)*hm3(i)-promro(i)
       g11=-(dim*ballmass(i)*hm4(i)+sumwh(i))
       g12=-sumwn(i)
       g2=indice(i)-nbaseline-deltahindex*(1.d0-2.d0/xplus)
       g21=-dnro(i)*sumwh(i)
       g22=(1.d0-dnro(i)*sumwn(i))
       ! Solving the matrix
       det=(g11*g22)-(g21*g12)
       detg1=(-g1*g22)-(-g2*g12)
       detg2=(-g2*g11)-(-g1*g21)
       deltah=detg1/det
       deltan=detg2/det
       if(l.le.liniNR) then
          deltah=0.d0
          deltan=0.d0
       endif
       h(i)=h(i)+deltah
       indice(i)=indice(i)+deltan
    enddo
!$omp end do

!$omp do schedule(static)
    do i=npini,npend
       call calculate_norm(indice(i),dim,pk(i),dpk(i))
    enddo
!$omp end do

!$omp end parallel

    
  END SUBROUTINE calculate_hNR_eq
