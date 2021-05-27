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
!            SPHYNX: calculate_avglogrho.f90          !
!                                                     !
! Calculates averaged log10(density) within 2h        !
!=====================================================!

  SUBROUTINE calculate_avglogrho

    USE parameters

    IMPLICIT NONE

    INTEGER i,iii,j,k

    if(flags)write(*,*) 'Average log density calculation'


    !Initialization
    roprom=0.d0      !Log Density

    !$omp parallel private(i,iii,k,j)
    !$omp do schedule(static)
    do i = npini,npend
       iii=i-npini+1
       do k=1,nvi(i)-1
          j=neighbors(iii,k)
          roprom(i)=roprom(i)+logprom(j)
       enddo
       roprom(i)=romprom(i)+logprom(i)   !Self-contribution
       roprom(i)=roprom(i)/dble(nvi(i))

    enddo
    !$omp end do
    !$omp end parallel

    RETURN
  END SUBROUTINE calculate_avglogrho
