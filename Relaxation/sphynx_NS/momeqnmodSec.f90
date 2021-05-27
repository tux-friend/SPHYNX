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
      DOUBLE PRECISION v1,v2,w1d,w2d,dw1,dw2,dter1,dter2
      DOUBLE PRECISION dpi,dpx,dpy,dpz,dpti,dpti0,dpxav,dpyav,dpzav
      DOUBLE PRECISION roti,fbalsi,divi,rotj,fbalsj,divj,fbals
      DOUBLE PRECISION vij1,vij2,vij3,vijrij,hij
      DOUBLE PRECISION muiji,tqij,dmy,dmy2,alfaij,betaij
      DOUBLE PRECISION rhoij,vijsignalu,dpxavu,dpyavu,dpzavu

      DOUBLE PRECISION kern11i,kern12i,kern13i,kern11j,kern12j,kern13j
      DOUBLE PRECISION kern21i,kern22i,kern23i,kern21j,kern22j,kern23j
      DOUBLE PRECISION kern31i,kern32i,kern33i,kern31j,kern32j,kern33j
      DOUBLE PRECISION termA1i,termA2i,termA3i,termA1j,termA2j,termA3j
      DOUBLE PRECISION dji1,dji2,dji3,wij,vijsignal,qiji,qijj,qijiu,qijju

      DOUBLE PRECISION,DIMENSION(n3)    :: aloc,vloc,gradploc,gradplocav
      DOUBLE PRECISION,DIMENSION(n)     :: promom,proene

      LOGICAL error
      
      call profile(0)
      if(flags)write(*,*) 'Moment & Energy equations (IAD0)'
      
      !Inicialization
      gradp=0.d0
      gradploc=0.d0
      gradplocav=0.d0
      energy=0.d0
      avisc=0.d0
      aviscu=0.d0
      muiji=0.d0
      maxvsignal=0.d0

      promom=pro
      proene=promom

      !grad-h terms
      if(gradh) then
         proene(:)=proene(:)/omega(:)
         promom(:)=promom(:)/omega(:)
      endif

!$omp parallel private(i,ii,ddmmax,k,j,jj,iii,d1,d2,d3,d02,d05,v1,v2,&
!$omp                  w1d,w2d,dw1,dw2,dter1,dter2,vij1,vij2,vij3,&
!$omp                  vijrij,divi,divj,roti,rotj,rhoij,vijsignalu,&
!$omp                  fbals,fbalsi,fbalsj,qiji,qijj,dpi,dpx,dpy,&
!$omp                  dpz,dpxav,dpyav,dpzav,tqij,kern11i,kern12i,&
!$omp                  kern13i,kern11j,kern12j,kern13j,kern21i,kern22i,&
!$omp                  kern23i,kern21j,kern22j,kern23j,kern31i,kern32i,&
!$omp                  kern33i,kern31j,kern32j,kern33j,termA1i,termA2i,&
!$omp                  termA3i,termA1j,termA2j,termA3j,dji1,dji2,dji3,&
!$omp                  jumpx,jumpy,jumpz,dmy,wij,vijsignal,qijiu,qijju,&
!$omp                  dpxavu,dpyavu,dpzavu,dmy2)


      !Reshaping
!$omp do schedule(static) 
      do i=1,n
         ii=1+dim*(i-1)
         aloc(ii)=a(i)
         aloc(ii+1)=a(i+n)
         aloc(ii+2)=a(i+n2)
         vloc(ii)=v(i)
         vloc(ii+1)=v(i+n)
         vloc(ii+2)=v(i+n2)
      enddo
!$omp end do


!$omp do schedule(static)              
         do i = npini,npend
            ii=1+dim*(i-1)
            iii=i-npini+1
            ddmmax = 0.d0

            do k=1,nvi(i)
               j=neighbors(iii,k)
               jj=1+dim*(j-1)

               call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

               if(i.ne.j)then
                  d1=aloc(ii)-aloc(jj)-jumpx
                  d2=aloc(ii+1)-aloc(jj+1)-jumpy
                  d3=aloc(ii+2)-aloc(jj+2)-jumpz
                  d02=d1*d1+d2*d2+d3*d3
                  d05=sqrt(d02)     
                  
                  v1=d05/h(i) 
                  v2=d05/h(j) 

                  call Wkernel(v1,indice(i),w1d,dw1)
                  call Wkernel(v2,indice(j),w2d,dw2)
                  
                  dter1=pk(i)*hm3(i)*w1d
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
                  
                  vij1=vloc(ii)-vloc(jj)
                  vij2=vloc(ii+1)-vloc(jj+1)
                  vij3=vloc(ii+2)-vloc(jj+2)
                  vijrij=vij1*d1+vij2*d2+vij3*d3

                  if(vijrij.lt.0.d0) then 
                     wij=vijrij/d05
                     rhoij=0.5d0*(promro(i)+promro(j))
                     vijsignal=c(i)+c(j)-3.d0*wij
                     vijsignalu=sqrt(abs(p(i)-p(j))/rhoij)
                     maxvsignal(i)=max(maxvsignal(i),vijsignal)

                     if(balsara) then
                        divi=abs(divv(i))
                        divj=abs(divv(j))
                        roti=abs(curlv(i))
                        rotj=abs(curlv(j))
                        fbalsi=divi/(divi+roti+5.d-4*c(i)/h(i))
                        fbalsj=divj/(divj+rotj+5.d-4*c(j)/h(j))
!                        if(fbalsi.le.0.05d0) fbalsi=0.05d0
!                        if(fbalsj.le.0.05d0) fbalsj=0.05d0
                     else
                        fbalsi=1.d0
                        fbalsj=1.d0
                     endif
                     qiji=-(alfa(i)+alfa(j))/4.d0*vijsignal*wij
                     qijiu=-alfau*vijsignalu*(u(i)-u(j))
                     qijj=qiji
                     qijju=qijiu
                     fbals=0.5*(fbalsi+fbalsj)
!                     fbals=2.d0*fbalsi*fbalsj/(fbalsi+fbalsj)
                     if(fbals.le.0.05d0) fbals=0.05d0
                     qiji=qiji*fbals
                     qijj=qijj*fbals
                  else 
                     qiji=0.d0  
                     qijj=0.d0
                     qijiu=0.d0  
                     qijju=0.d0
                  endif

                  dpx=promom(i)*termA1i+promom(j)*termA1j
                  dpy=promom(i)*termA2i+promom(j)*termA2j
                  dpz=promom(i)*termA3i+promom(j)*termA3j
                  dmy=vol(i)*masa(j)/masa(i)*qiji
                  dmy2=vol(j)*qijj
                  dpxav=0.5d0*(dmy*termA1i+dmy2*termA1j)
                  dpyav=0.5d0*(dmy*termA2i+dmy2*termA2j)
                  dpzav=0.5d0*(dmy*termA3i+dmy2*termA3j)
                  dmy=vol(i)*masa(j)/masa(i)*qijiu
                  dmy2=vol(j)*qijju
                  dpxavu=0.5d0*(dmy*termA1i+dmy2*termA1j)
                  dpyavu=0.5d0*(dmy*termA2i+dmy2*termA2j)
                  dpzavu=0.5d0*(dmy*termA3i+dmy2*termA3j)
                  
                  dpi=proene(i)*2.d0
                  
                  tqij=vij1*termA1i+vij2*termA2i+vij3*termA3i
                  
                  dmy=xmass(i)*xmass(j)/masa(i) 
                  gradploc(ii)=gradploc(ii)+dmy*dpx+dpxav
                  gradploc(ii+1)=gradploc(ii+1)+dmy*dpy+dpyav
                  gradploc(ii+2)=gradploc(ii+2)+dmy*dpz+dpzav
                  energy(i)=energy(i)+dmy*dpi*tqij
                  avisc(i)=avisc(i)+dpxav*vij1+dpyav*vij2+dpzav*vij3
                  aviscu(i)=aviscu(i)+dpxavu*vij1+dpyavu*vij2+dpzavu*vij3
               endif
               
            enddo
            
         enddo
!$omp end do

      !Reshaping back
!$omp do schedule(static)
      do i=1,n
         ii=1+dim*(i-1)
         gradp(i)=gradploc(ii)
         gradp(i+n)=gradploc(ii+1)
         gradp(i+n2)=gradploc(ii+2)
         avisc(i)=max(0.d0,avisc(i))
!         aviscu(i)=max(0.d0,aviscu(i))
      enddo
!$omp end do

!$omp end parallel 

      call profile(1)

      RETURN
    END SUBROUTINE momeqn
