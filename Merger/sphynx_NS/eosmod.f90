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
!                  SPHYNX: eosmod.f90                 !
!                                                     !
! Calls different EOS.                                !
!=====================================================!

      SUBROUTINE eostot

      USE parameters
      USE eosmodule

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(nmax):: ueos

      INTEGER i,keytemp,keyerr,m,mitermax,itov,o
      DOUBLE PRECISION xrho,xye,xtemp,xtemp2
      DOUBLE PRECISION xenr,xprs,xent,xcs2,xdedt,xmunu
      DOUBLE PRECISION xdpderho,xdpdrhoe
      DOUBLE PRECISION iye,yel,yeg,dmutot
      DOUBLE PRECISION xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
      DOUBLE PRECISION xxa,xxh,xxn,xxp,Tl,Tg

      call profile(0)

      if(flags)print *,'Calculating EOS'

      p=0.d0
      c=0.d0
      pro=0.d0
      dpdt=0.d0
      cv=0.d0
      dudv=0.d0
      ueos=0.d0

      if (eos.eq.1) then

         call eosid_t(ueos)
         !input: temp,promro,mui
         !This populates p,u,c,cv,dpdt,dpdr

      else if (eos.eq.2) then

         call eosid_u
         !input: promro,u,gamma
         !This populates p,temp,c,cv,dpdt
         if(l.le.liniNR.and.KH) then
            p(:)=p0
            u(:)=p0/promro(:)/(gamma-1.d0)
            c(:)=sqrt((gamma-1.d0)*u(:))
         endif

      else if (eos.eq.3) then
         p(:)=kpol*promro(:)**gammapol
         c(:)=sqrt(gammapol*p(:)/promro(:))

      else if (eos.eq.4) then
         if(l.eq.1)then
             keytemp=1
             keyerr=0
             xtemp=0.5d0
             do i=1,nmax
                o=0
                yel=eos_yemin		
	            yeg=eos_yemax
	            yelec(i)=(yel+yeg)*0.5d0
	            call nuc_eos_full(promro(i),xtemp,yelec(i),xenr,xprs,xent,xcs2,xdedt,&
		             xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
		             xmuhat,keytemp,keyerr,precision)
	            dmutot=xmu_p+xmu_e-xmu_n  
	            do while(abs(dmutot).ge.precision.and.(o.le.50))
	               if(dmutot.le.precision)then
	                  yel=yelec(i)
	                  yelec(i)=(yelec(i)+yeg)*0.5d0
                   else
                      yeg=yelec(i)
			          yelec(i)=(yelec(i)+yel)*0.5d0
		          endif
		          call nuc_eos_full(promro(i),xtemp,yelec(i),xenr,xprs,xent,xcs2,xdedt,&
			           xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,&
			           xmu_p,xmuhat,keytemp,keyerr,precision)
		          dmutot=xmu_p+xmu_e-xmu_n
		          o=o+1
	           enddo		   
			   call nuc_eos_full(promro(i),temp(i),yelec(i),xenr,p(i),xent,xcs2,xdedt,&
		              xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,&
			      xmu_p,xmuhat,keytemp,keyerr,precision)
			   c(i)=sqrt(xcs2)
			   u(i)=xenr
			enddo
         else
			do i=1,nmax
			   keytemp=1
			   keyerr=0
			   xtemp=0.5d0
	           o=0
                  yel=eos_yemin		
	           yeg=eos_yemax
	           yelec(i)=(yel+yeg)*0.5d0
	           call nuc_eos_full(promro(i),xtemp,yelec(i),xenr,xprs,xent,xcs2,xdedt,&
		            xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
		            xmuhat,keytemp,keyerr,precision)
	           dmutot=xmu_p+xmu_e-xmu_n  
	           do while(abs(dmutot).ge.precision.and.(o.le.50))
		          if(dmutot.le.precision)then
			         yel=yelec(i)
			         yelec(i)=(yelec(i)+yeg)*0.5d0
		          else	
			         yeg=yelec(i)
			         yelec(i)=(yelec(i)+yel)*0.5d0
		          endif
		          call nuc_eos_full(promro(i),xtemp,yelec(i),xenr,xprs,xent,xcs2,xdedt,&
			           xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,&
			           xmu_p,xmuhat,keytemp,keyerr,precision)
		          dmutot=xmu_p+xmu_e-xmu_n
		          o=o+1
	           enddo
	           keytemp=0
	           !print *,'particle:',i,promro(i),yelec(i),u(i),temp(i)
			   call nuc_eos_full(promro(i),temp(i),yelec(i),u(i),p(i),xent,xcs2,xdedt,&
			      xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
		          xmuhat,keytemp,keyerr,precision)
		       c(i)=sqrt(xcs2)
			enddo
		 endif
      else if (eos.eq.5) then
         do i=1,nmax
            if(promro(i).lt.rhoth0)then
                p(i)=kpol0*promro(i)**gammapol0
                c(i)=sqrt(gammapol0*p(i)/promro(i))
            else if((promro(i).lt.rhoth1).and.(promro(i).ge.rhoth0))then
                p(i)=kpol1*promro(i)**gammapol1
                c(i)=sqrt(gammapol1*p(i)/promro(i))
            else if((promro(i).lt.rhoth2).and.(promro(i).ge.rhoth1))then
                p(i)=kpol2*promro(i)**gammapol2
                c(i)=sqrt(gammapol2*p(i)/promro(i))
            else
                p(i)=kpol3*promro(i)**gammapol3
                c(i)=sqrt(gammapol3*p(i)/promro(i))
            endif
         enddo
      else

         write(*,*) 'Not defined EOS!'
         stop

      endif

      cmax=maxval(c)
      !if(l.eq.iterini.and.inienergy)u=ueos

      if(std_VE) then
         pro(:)=p(:)!/(sumkx(:))**2
      else
         pro(:)=p(:)/(sumkx(:))
      endif

      call profile(1)

      RETURN
      END SUBROUTINE eostot
