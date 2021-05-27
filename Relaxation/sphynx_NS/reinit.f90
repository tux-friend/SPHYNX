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
!                  SPHYNX: reinit.f90                 !
!                                                     !
! Re-initialization of several arrays.                !
!=====================================================!

      subroutine reinit()
        
      use parameters

      implicit none

      !For integration eq 1 this is directly done in actualizamod.
      u0=u
      temp0=temp
      s0=s
      ye0=ye
!      xn0=xn
      nuedot=0.d0
      
      return
      end subroutine reinit
