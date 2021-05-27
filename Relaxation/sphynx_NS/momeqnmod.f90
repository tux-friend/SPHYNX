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
!                SPHYNX: momeqnmod.f90                !
!                                                     !
! Calculates momentum and energy equations including  !
! IAD terms as explained in                           !
! Garcia-Senz, Cabezon, Escartin, A&A, 538, A9 (2012) !
!                                                     !
!=====================================================!

      SUBROUTINE momeqn

      USE parameters

      IMPLICIT NONE

      INTEGER i,j,k,ii,jj,iii

      DOUBLE PRECISION d1,d2,d3,d05,d02,jumpx,jumpy,jumpz
      DOUBLE PRECISION v1,v2,w1d,w2d,dter1,dter2,norm
      DOUBLE PRECISION dpi,dpx,dpy,dpz,dpti,dpti0,dpxav,dpyav,dpzav
      DOUBLE PRECISION roti,fbalsi,divi,rotj,fbalsj,divj,fbals
      DOUBLE PRECISION vij1,vij2,vij3,vijrij,hij
      DOUBLE PRECISION muiji,tqij,dmy,dmy1,dmy2,alfaij,betaij
      DOUBLE PRECISION rhoij,vijsignalu,dpxavu,dpyavu,dpzavu
      DOUBLE PRECISION uijx,uijy,uijz,At_ab,sigma_ab

      DOUBLE PRECISION kern11i,kern12i,kern13i,kern11j,kern12j,kern13j
      DOUBLE PRECISION kern21i,kern22i,kern23i,kern21j,kern22j,kern23j
      DOUBLE PRECISION kern31i,kern32i,kern33i,kern31j,kern32j,kern33j
      DOUBLE PRECISION termA1i,termA2i,termA3i,termA1j,termA2j,termA3j
      DOUBLE PRECISION dji1,dji2,dji3,wij,vijsignal,qiji,qijj,qijiu,qijju

      DOUBLE PRECISION,DIMENSION(n)     :: promom,proene

      LOGICAL error

      call profile(0)
      if(flags)write(*,*) 'Moment & Energy equations (IAD0)'

      !Inicialization
      gradp=0.d0
      energy=0.d0
      avisc=0.d0
      aviscu=0.d0
      muiji=0.d0
      maxvsignal=0.d0
      mark_ramp=0.d0

      promom=pro
      proene=promom

!       if(std_VE)then
!          if(.not.crossedXmass) then
!              proene(:)=proene(:)/sumkx(:)**2
!              promom(:)=promom(:)/sumkx(:)**2
!          endif
!       else
! #ifdef EQMASS
!          if(crossedXmass) then
!             proene(:)=proene(:)/masspart**2
!             promom(:)=promom(:)/masspart**2
!          else
!             proene(:)=proene(:)/masspart**2*xmass(:)**2
!             promom(:)=promom(:)/masspart**2*xmass(:)**2
!          endif
! #else
!          if(crossedXmass) then
!             proene(:)=proene(:)/mass(:)**2
!             promom(:)=promom(:)/mass(:)**2
!          else
!             proene(:)=proene(:)/mass(:)**2*xmass(:)**2
!             promom(:)=promom(:)/mass(:)**2*xmass(:)**2
!          endif
! #endif
!       endif

      if(.not.std_VE) then
         proene(:)=proene(:)/masspart**2
         promom(:)=promom(:)/masspart**2
      endif


      !grad-h terms
      if(gradh) then
          proene(:)=proene(:)/omega(:)
          promom(:)=promom(:)/omega(:)
      endif

      !$omp parallel private(i,ii,ddmmax,k,j,jj,iii,d1,d2,d3,d02,d05,v1,v2,&
      !$omp                  w1d,w2d,dter1,dter2,norm,vij1,vij2,vij3,&
      !$omp                  vijrij,divi,divj,roti,rotj,rhoij,vijsignalu,&
      !$omp                  fbals,fbalsi,fbalsj,qiji,qijj,dpi,dpx,dpy,&
      !$omp                  dpz,dpxav,dpyav,dpzav,tqij,kern11i,kern12i,&
      !$omp                  kern13i,kern11j,kern12j,kern13j,kern21i,kern22i,&
      !$omp                  kern23i,kern21j,kern22j,kern23j,kern31i,kern32i,&
      !$omp                  kern33i,kern31j,kern32j,kern33j,termA1i,termA2i,&
      !$omp                  termA3i,termA1j,termA2j,termA3j,dji1,dji2,dji3,&
      !$omp                  jumpx,jumpy,jumpz,dmy,wij,vijsignal,qijiu,qijju,&
      !$omp                  dpxavu,dpyavu,dpzavu,dmy1,dmy2,uijx,uijy,uijz,&
      !$omp                  At_ab,sigma_ab)


      !$omp do schedule(static)
      do i = npini,npend
         ii=1+dim*(i-1)
         iii=i-npini+1
         ddmmax = 0.d0
         norm=pk(i)*hm3(i)

         do k=1,nvi(i)-1
            j=neighbors(iii,k)
            jj=1+dim*(j-1)

            call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

            d1=a(ii)-a(jj)-jumpx
            d2=a(ii+1)-a(jj+1)-jumpy
            d3=a(ii+2)-a(jj+2)-jumpz
            d02=d1*d1+d2*d2+d3*d3
            d05=sqrt(d02)

            uijx=d1/d05
            uijy=d2/d05
            uijz=d3/d05

            v1=d05/h(i)
            v2=d05/h(j)

            call Wkernel_noderiv(v1,indice(i),w1d)
            call Wkernel_noderiv(v2,indice(j),w2d)

            dter1=norm*w1d
            dter2=pk(j)*hm3(j)*w2d

            dji1=-d1
            dji2=-d2
            dji3=-d3

            kern11i=c11(i)*dji1
            kern12i=c12(i)*dji2
            kern13i=c13(i)*dji3
            kern11j=c11(j)*dji1
            kern12j=c12(j)*dji2
            kern13j=c13(j)*dji3

            kern21i=c12(i)*dji1
            kern22i=c22(i)*dji2
            kern23i=c23(i)*dji3
            kern21j=c12(j)*dji1
            kern22j=c22(j)*dji2
            kern23j=c23(j)*dji3

            kern31i=c13(i)*dji1
            kern32i=c23(i)*dji2
            kern33i=c33(i)*dji3
            kern31j=c13(j)*dji1
            kern32j=c23(j)*dji2
            kern33j=c33(j)*dji3

            termA1i=(kern11i+kern12i+kern13i)*dter1
            termA2i=(kern21i+kern22i+kern23i)*dter1
            termA3i=(kern31i+kern32i+kern33i)*dter1
            termA1j=(kern11j+kern12j+kern13j)*dter2
            termA2j=(kern21j+kern22j+kern23j)*dter2
            termA3j=(kern31j+kern32j+kern33j)*dter2


!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
!                                                                            !
!                             ARTIFICIAL VISCOSITY                           !
!                                                                            !
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
!                                                                            !
!     Monaghan (1992) THE VISCOSITY VANISHES WHEN vijrij>0, WHICH IS         !
!     THE SPH EQUIVALENT OF THE CONDITION V.v>0. THE EXPRESSION qij          !
!     CONTAINS A TERM THAT IS LINEAR IN THE VELOCITY DIFFERENCES, WHICH      !
!     PRODUCES A SHEAR AND BULK VISCOSITY. THE QUADRATIC TERM IS NECESSARY   !
!     TO HANDLE HIGH MACH NUMBER SHOCKS, AND IS ROUGHLY EQUIVALENT TO THE    !
!     VON NEUMANN-RICHTMYER VISCOSITY USED IN FINITE-DIFFERENCE METHODS ===> !
!     ===> it is Galilean invariant - it vanishes for rigid rotation - it    !
!     conserves total linear and angular momenta.                            !
!                                                                            !
!----------------------------------------------------------------------------!

            vij1=v(ii)-v(jj)
            vij2=v(ii+1)-v(jj+1)
            vij3=v(ii+2)-v(jj+2)
            vijrij=vij1*d1+vij2*d2+vij3*d3

            if(vijrij.lt.0.d0) then
               wij=vijrij/d05
               rhoij=0.5d0*(promro(i)+promro(j))
               !vijsignal=(alfa(i)+alfa(j))/4.d0*(c(i)+c(j))-betaAV*wij
               vijsignal=c(i)+c(j)-3.d0*wij
               vijsignalu=sqrt(abs(p(i)-p(j))/rhoij)
               maxvsignal(i)=max(maxvsignal(i),vijsignal)
               vijsignal=(alfa(i)+alfa(j))/4.d0*(c(i)+c(j))-betaAV*wij

               if(balsara) then
                  divi=abs(divv(i))
                  divj=abs(divv(j))
                  roti=abs(curlv(i))
                  rotj=abs(curlv(j))
                  fbalsi=divi/(divi+roti+5.d-4*c(i)/h(i))
                  fbalsj=divj/(divj+rotj+5.d-4*c(j)/h(j))
                  ! if(fbalsi.le.0.05d0) fbalsi=0.05d0
                  ! if(fbalsj.le.0.05d0) fbalsj=0.05d0
               else
                  fbalsi=1.d0
                  fbalsj=1.d0
               endif
               qiji=-vijsignal*wij
               !qiji=-(alfa(i)+alfa(j))/4.d0*vijsignal*wij
               qijiu=alfau*vijsignalu*(u(i)-u(j))
               qijj=qiji
               qijju=qijiu
               fbals=0.5d0*(fbalsi+fbalsj)
                !  fbals=2.d0*fbalsi*fbalsj/(fbalsi+fbalsj)
               if(fbals.le.0.05d0) fbals=0.05d0
               qiji=qiji*fbals
               qijj=qijj*fbals
            else
               qiji=0.d0
               qijj=0.d0
               qijiu=0.d0
               qijju=0.d0
            endif

            At_ab=(abs(promro(i)-promro(j)))/(promro(i)+promro(j))
#ifdef EQMASS
            if(std_VE) then
               if(At_ab.lt.Atmin) then
                  dmy1=xmass(i)*xmass(j)/masspart/sumkx(i)**2
                  dmy2=xmass(i)*xmass(j)/masspart/sumkx(j)**2
               else if(At_ab.gt.Atmax) then
                  dmy1=xmass(i)*xmass(j)/masspart/sumkx(i)/sumkx(j)
                  dmy2=dmy1
                  mark_ramp(i)=mark_ramp(i)+1.d0
               else
                  sigma_ab=ramp*(At_ab-Atmin)
                  dmy1=xmass(i)*xmass(j)/masspart/sumkx(i)**(2.d0-sigma_ab)/sumkx(j)**sigma_ab
                  dmy2=xmass(i)*xmass(j)/masspart/sumkx(j)**(2.d0-sigma_ab)/sumkx(i)**sigma_ab
                  mark_ramp(i)=mark_ramp(i)+sigma_ab
               endif

            else
               if(At_ab.lt.Atmin) then
                  dmy1=masspart*xmass(i)**2
                  dmy2=masspart*xmass(j)**2
               else if(At_ab.gt.Atmax) then
                  dmy1=masspart*xmass(i)*xmass(j)
                  dmy2=dmy1
                  mark_ramp(i)=mark_ramp(i)+1.d0
               else
                  sigma_ab=ramp*(At_ab-Atmin)
                  dmy1=masspart*xmass(i)**(2.d0-sigma_ab)*xmass(j)**sigma_ab
                  dmy2=masspart*xmass(j)**(2.d0-sigma_ab)*xmass(i)**sigma_ab
                  mark_ramp(i)=mark_ramp(i)+sigma_ab
               endif

            endif

#else
            if(std_VE) then
              if(At_ab.lt.Atmin) then
                dmy1=xmass(i)*xmass(j)/mass(i)/sumkx(i)**2
                dmy2=xmass(i)*xmass(j)/mass(i)/sumkx(j)**2
              else if(At_ab.gt.Atmax) then
                dmy1=xmass(i)*xmass(j)/mass(i)/sumkx(i)/sumkx(j)
                dmy2=dmy1
                mark_ramp(i)=mark_ramp(i)+1.d0
              else
                sigma_ab=ramp*(At_ab-Atmin)
                dmy1=xmass(i)*xmass(j)/mass(i)/sumkx(i)**(2.d0-sigma_ab)/sumkx(j)**sigma_ab
                dmy2=xmass(i)*xmass(j)/mass(i)/sumkx(j)**(2.d0-sigma_ab)/sumkx(i)**sigma_ab
                mark_ramp(i)=mark_ramp(i)+sigma_ab
              endif

            else
              if(At_ab.lt.Atmin) then
                dmy1=mass(j)*xmass(i)**2
                dmy2=mass(j)*xmass(j)**2
              else if(At_ab.gt.Atmax) then
                dmy1=mass(j)*xmass(i)*xmass(j)
                dmy2=dmy1
                mark_ramp(i)=mark_ramp(i)+1.d0
              else
                sigma_ab=ramp*(At_ab-Atmin)
                dmy1=mass(j)*xmass(i)**(2.d0-sigma_ab)*xmass(j)**sigma_ab
                dmy2=mass(j)*xmass(j)**(2.d0-sigma_ab)*xmass(i)**sigma_ab
                mark_ramp(i)=mark_ramp(i)+sigma_ab
              endif

           endif
#endif
            dpx=promom(i)*dmy1*termA1i+promom(j)*dmy2*termA1j
            dpy=promom(i)*dmy1*termA2i+promom(j)*dmy2*termA2j
            dpz=promom(i)*dmy1*termA3i+promom(j)*dmy2*termA3j

!             dpx=promom(i)*termA1i+promom(j)*termA1j
!             dpy=promom(i)*termA2i+promom(j)*termA2j
!             dpz=promom(i)*termA3i+promom(j)*termA3j
#ifdef EQMASS
            dmy=vol(i)*qiji
#else
            dmy=vol(i)*mass(j)/mass(i)*qiji
#endif
            dmy2=vol(j)*qijj
            dpxav=0.5d0*(dmy*termA1i+dmy2*termA1j)
            dpyav=0.5d0*(dmy*termA2i+dmy2*termA2j)
            dpzav=0.5d0*(dmy*termA3i+dmy2*termA3j)
#ifdef EQMASS
            dmy=vol(i)*qijiu
#else
            dmy=vol(i)*mass(j)/mass(i)*qijiu
#endif
            dmy2=vol(j)*qijju
            dpxavu=0.5d0*(dmy*termA1i+dmy2*termA1j)
            dpyavu=0.5d0*(dmy*termA2i+dmy2*termA2j)
            dpzavu=0.5d0*(dmy*termA3i+dmy2*termA3j)

            dpi=proene(i)*2.d0

            tqij=vij1*termA1i+vij2*termA2i+vij3*termA3i

! #ifdef EQMASS
!             if(std_VE) then
!                if(crossedXmass) then
!                   dmy=xmass(i)*xmass(j)/masspart/sumkx(i)/sumkx(j)
!                else
!                   dmy=xmass(i)*xmass(j)/masspart
!                endif
!             else
!                if(crossedXmass) then
!                   dmy=masspart*xmass(i)*xmass(j)
!                else
!                   dmy=masspart
!                endif
!             endif
! #else
!             if(std_VE) then
!                if(crossedXmass) then
!                   dmy=xmass(i)*xmass(j)/mass(i)/sumkx(i)/sumkx(j)
!                else
!                   dmy=xmass(i)*xmass(j)/mass(i)
!                endif
!             else
!                if(crossedXmass) then
!                   dmy=mass(j)*xmass(i)*xmass(j)
!                else
!                   dmy=mass(j)
!                endif
!             endif
! #endif
!            gradp(ii)=gradp(ii)+dmy*dpx+dpxav
!            gradp(ii+1)=gradp(ii+1)+dmy*dpy+dpyav
!            gradp(ii+2)=gradp(ii+2)+dmy*dpz+dpzav
            gradp(ii)=gradp(ii)+dpx+dpxav
            gradp(ii+1)=gradp(ii+1)+dpy+dpyav
            gradp(ii+2)=gradp(ii+2)+dpz+dpzav
            !energy(i)=energy(i)+dmy*dpi*tqij
            energy(i)=energy(i)+dmy1*dpi*tqij
            avisc(i)=avisc(i)+dpxav*vij1+dpyav*vij2+dpzav*vij3
            if(isnan(avisc(i))) then
              print *, i,j,nvi(i)
              print *, dpxav,dpyav,dpzav
              print *, vij1,vij2,vij3
              print *, qiji,qijj
              print *, vol(i),vol(j)
              print *, sigma_ab,mark_ramp(i)
              print *, promro(i),promro(j)
              print *, h(i),h(j)
              print *, c(i),c(j)
              print *, p(i),p(j)
              print *, alfa(i),alfa(j)
              stop
            endif
            !aviscu(i)=aviscu(i)+dpxavu*vij1+dpyavu*vij2+dpzavu*vij3
            aviscu(i)=aviscu(i)+dpxavu*uijx+dpyavu*uijy+dpzavu*uijz
         enddo
      enddo
      !$omp end do
      !$omp end parallel

      mark_ramp(:)=mark_ramp(:)/dble(nvi(:))
      avisc(:)=max(0.d0,avisc(:))
! #ifdef EQMASS
!       if(std_VE) then
!          gradp(:)=gradp(:)/masspart
!          energy(:)=energy(:)/masspart
!       else
!          gradp(:)=gradp(:)*masspart
!          energy(:)=energy(:)*masspart
!       endif
! #endif


      call profile(1)

      RETURN
    END SUBROUTINE momeqn
