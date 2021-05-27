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
!                SPHYNX: estabilmod.f90               !
!                                                     !
! Calculates central density.                         !
!=====================================================!

   SUBROUTINE estabilmod()
       
     USE parameters

     IMPLICIT NONE
     INTEGER i,j,k,ncount,idmy,kminb
     INTEGER,DIMENSION(n) :: indx2
     DOUBLE PRECISION rrrmax,ax,ay,az,rrr,dmy,rhocloud,uambient,masacloud,L1unit,L1deltar
     DOUBLE PRECISION xcmb,ycmb,zcmb,dmin,masab

     !find promedium central density and max radius of star (with 50 part.)
     call indexx(n,promro,indx2)
     rocentral=0.d0
     rrrmax=0.d0
     refax=0.d0
     refay=0.d0
     refaz=0.d0
     do i=n,n-49,-1  !Tiene que ser 49 para que hayan 50 particulas.
        j=indx(i)
        k=indx2(i)
        rrrmax=rrrmax+radius(j)
        rocentral=rocentral+promro(k)
        refax=refax+a(k)-despl
        refay=refay+a(k+n)-despl
        refaz=refaz+a(k+n2)-despl
     enddo
     rocentral=rocentral/50.d0
     rrrmax=rrrmax/50.d0
     refax=refax/50.d0
     refay=refay/50.d0
     refaz=refaz/50.d0

!     Center of mass Blob
      xcmb=0.d0
      ycmb=0.d0
      zcmb=0.d0 
        masab=0.d0
      do i=1,n
        if(mark(i).eq.1.d0) then 
          xcmb=xcmb+a(i)*masa(i)
          ycmb=ycmb+a(i+n)*masa(i)
          zcmb=zcmb+a(i+n2)*masa(i)
          masab=masab+masa(i)
        endif
      enddo
       xcmb=xcmb/masab
       ycmb=ycmb/masab
       zcmb=zcmb/masab
       dmin=1.d80
!       find closest particle to center of mass blob
       do i=1,n
        if(mark(i).eq.1.d0) then 
        dmy=sqrt((a(i)-xcmb)**2+(a(i+n)-ycmb)**2+(a(i+n2)-zcmb)**2)
          if(dmy.le.dmin) then 
            dmin=dmy
           kminb=i
        endif
       endif
       enddo
       

!     Calculation of L1unit,L1deltar. Only for particles at distance D from particle 
!   at cdm bubble

      L1unit=0.d0
      L1deltar=0.d0
      ncount=0

      do i=1,n
        dmy=sqrt((a(i)-a(kminb))**2+(a(i+n)-a(kminb+n))**2+(a(i+n2)-a(kminb+n2))**2)
       if(dmy.le.0.075) then 
        dmy=sqrt(checkdeltax(i)**2+checkdeltay(i)**2+checkdeltaz(i)**2)
        L1unit=L1unit+abs(checkInorm(i)-1.d0)
        L1deltar=L1deltar+dmy/h(i)
        ncount=ncount+1
      endif
      enddo
        L1unit=L1unit/ncount
        L1deltar=L1deltar/ncount

     
!     Calculation of stripped mass of the blob

         rhocloud=10.d0
         uambient=3.d0/2.d0
         masacloud=0.0d0
      do i=1,n
        if(promro(i).ge.0.64*rhocloud .and. u(i).le.0.9*uambient) then
           masacloud=masacloud+masa(i)
        endif
      enddo
         
     if(id.eq.0)then
        open(11,file='estabil.d',access='append')
        rrr=sqrt(refax**2+refay**2+refaz**2)
        write(11,formatin) tt,rocentral,masacloud/masacloudinic,L1unit,L1deltar,real(ncount),xcmb-despl,ycmb-despl,zcmb-despl,promro(kminb),p(kminb)
        close(11)
     endif

     if(rocentral.lt.1.d12)then
        refrad=0.d0
        refax=0.d0
        refay=0.d0
        refaz=0.d0
     else
        do i=1,n
           ax=a(i)-despl-refax
           ay=a(i+n)-despl-refay
           az=a(i+n2)-despl-refaz
           radius(i)=sqrt(ax*ax+ay*ay+az*az)
        enddo
     endif

     call indexx(n,radius,indx)  !Indx is sorted by radius measured from
                                 !0,0 when rho<1e12 and from where density
                                 !is higher if rho>1e12.

     RETURN
   END SUBROUTINE
