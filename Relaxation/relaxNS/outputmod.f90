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
!                SPHYNX: outputmod.f90                !
!                                                     !
! Output of data.                                     !
!=====================================================!

    SUBROUTINE output

      USE parameters

      IMPLICIT NONE
      INTEGER i,j,k,ii
      DOUBLE PRECISION macum
      DOUBLE PRECISION,DIMENSION(dim)::ps
      DOUBLE PRECISION,DIMENSION(nmax)::vrad,gpr,fgr
      CHARACTER*17 nomfs
      CHARACTER*9 nomfs1
      CHARACTER*3 prefix1
      CHARACTER*6 sufix
      CHARACTER*14 arch

      DOUBLE PRECISION masscloud
      DOUBLE PRECISION,PARAMETER:: rhocloud=10.d0
      DOUBLE PRECISION,PARAMETER:: uambient=1.5d0
      DOUBLE PRECISION,PARAMETER:: tkh=0.0937d0


      call profile(0)

      !user-defined conditions for writing
      if(mod(l,iodelay).eq.0.or.(l.eq.iterini.and.outini))escribe=.true.
      if(checkdens)escribe=.true.
      if(l.eq.86)escribe=.true.

      if(flags.and.escribe)write(*,*)'Writing data'

      if(escribe) then
         open(10,file='outputtimes.d',position='append')
         write(10,'(1x,i6,1x,a17,2(1x,es17.10))') l,'                 ',tt
         close(10)
      endif

      if(escribe) then
         escribe=.false.
         prefix1='s1.'

         call nomfils(nomfs1,l,prefix1)
         nomfs='data/'//nomfs1
         write(*,*) 'Output file: ',nomfs1
         open(10,file=nomfs)!,form='unformatted')

!--------------------  User can change this  --------------------
         macum=0.d0
         do i=1,n1
            k=i
            ii=1+dim*(k-1)
            ps(1)=a(ii)-despl
            ps(2)=a(ii+1)-despl
            ps(3)=a(ii+2)-despl
            vrad(k)=(ps(1)*v(ii)+ps(2)*v(ii+1)+ps(3)*v(ii+2))/radius(k)
            gpr(k)=(ps(1)*gradp(ii)+ps(2)*gradp(ii+1)+ps(3)*gradp(ii+2))/radius(k)
            fgr(k)=g*(ps(1)*f(ii)+ps(2)*f(ii+1)+ps(3)*f(ii+2))/radius(k)
#ifdef EQMASS
            macum=macum+masspart
#else
            macum=macum+mass(k)
#endif
            write(10,formatout) k,(ps(j),j=1,dim),h(k),u(k),promro(k),&
                 &    v(ii),v(ii+1),v(ii+2),radius(k),yelec(k),p(k),&
                 &    divv(k),omega(k),f(ii),f(ii+1),f(ii+2),gradp(ii),gradp(ii+1),&
                 &    gradp(ii+2),temp(k),avisc(k)*0.5d0,energy(k)*.5d0,&
                 &    ballmass(k),xmass(k),vrad(k),dble(nvi(k)),alfa(k),&
                 &    mark_ramp(k),maxvsignal(k),c(k),sumkx(k),checkInorm(k),&
                 &    fgr(k),gpr(k)
         enddo
!----------------------------------------------------------------

         close(10)

      endif



      call profile(1)

      RETURN
    END SUBROUTINE output
