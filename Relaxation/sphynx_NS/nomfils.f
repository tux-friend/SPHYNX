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
!                 SPHYNX: nomfils.f90                 !
!                                                     !
! Creates a name for a file output depending on an    !
! integer.                                            !
!=====================================================!

      SUBROUTINE nomfils(nomfs,i,prefix)
      implicit none
      integer,INTENT(in)::i
      integer centmillar,decmillar,millar,cent,dec,uni,rest
      character,INTENT(out)::nomfs*9 
      character,INTENT(in)::prefix*3
      character sufix*6
      centmillar=i/100000+48
      rest=mod(i,100000)
      decmillar=rest/10000+48
      rest=mod(i,10000)
      millar=rest/1000+48
      rest=mod(i,1000)
      cent=rest/100+48
      rest=mod(i,100)
      dec=rest/10+48
      uni=mod(i,10)+48

      sufix=char(centmillar)//char(decmillar)//char(millar)//char(cent)
     *     //char(dec)//char(uni)
      nomfs=prefix//sufix
      RETURN
      END
