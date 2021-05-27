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
     INTEGER i,j,k,ii,jj,ncount,kminb
     INTEGER,DIMENSION(n) :: indx2
     DOUBLE PRECISION rrrmax,ax,ay,az,rrr
     DOUBLE PRECISION xcmb,ycmb,zcmb,masab,dmin,dmy,L1unit,L1deltar

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
        ii=1+dim*(k-1)
        rrrmax=rrrmax+radius(j)
        rocentral=rocentral+promro(k)
        refax=refax+a(ii)-despl
        refay=refay+a(ii+1)-despl
        refaz=refaz+a(ii+2)-despl
     enddo
     rocentral=rocentral/50.d0
     rrrmax=rrrmax/50.d0
     refax=refax/50.d0
     refay=refay/50.d0
     refaz=refaz/50.d0

     !Center of mass Blob
     xcmb=0.d0
     ycmb=0.d0
     zcmb=0.d0
     masab=0.d0
     do i=1,n
       if(mark(i)) then
         ii=1+dim*(i-1)
         xcmb=xcmb+a(ii)
         ycmb=ycmb+a(ii+1)
         zcmb=zcmb+a(ii+2)
         masab=masab+1.d0
       endif
     enddo
     xcmb=xcmb/masab
     ycmb=ycmb/masab
     zcmb=zcmb/masab
     dmin=1.d80
     kminb=0

     !find closest particle to center of mass blob
     do i=1,n
       if(mark(i)) then
         ii=1+dim*(i-1)
         dmy=sqrt((a(ii)-xcmb)**2+(a(ii+1)-ycmb)**2+(a(ii+2)-zcmb)**2)
         if(dmy.le.dmin) then
           dmin=dmy
           kminb=i
         endif
       endif
     enddo
     if(kminb.lt.1.or.kminb.gt.nmax) print *,"Wrong Kminb in estabilmod.f90",kminb
     !Calculation of L1unit,L1deltar. Only for particles at distance D from particle
     !at cdm bubble
     L1unit=0.d0
     L1deltar=0.d0
     ncount=0
     jj=1+dim*(kminb-1)

     do i=1,n
       ii=1+dim*(i-1)
       dmy=sqrt((a(ii)-a(jj))**2+(a(ii+1)-a(jj+1))**2+(a(ii+2)-a(jj+2))**2)
       if(dmy.le.0.075d0) then
         dmy=sqrt(checkdeltax(i)**2+checkdeltay(i)**2+checkdeltaz(i)**2)
         L1unit=L1unit+abs(checkInorm(i)-1.d0)
         L1deltar=L1deltar+dmy/h(i)
         ncount=ncount+1
       endif
     enddo
     L1unit=L1unit/dble(ncount)
     L1deltar=L1deltar/dble(ncount)


     if(id.eq.0)then
        open(11,file='estabil.d',access='append')
        rrr=sqrt(refax**2+refay**2+refaz**2)
        write(11,formatin) tt,rrrmax,rocentral,L1unit,L1deltar,dble(ncount),&
              & xcmb-despl,ycmb-despl,zcmb-despl,promro(kminb),p(kminb)
        close(11)
     endif

     if(rocentral.lt.1.d12)then
        refrad=0.d0
        refax=0.d0
        refay=0.d0
        refaz=0.d0
     else
        do i=1,n
           ii=1+dim*(i-1)
           ax=a(ii)-despl-refax
           ay=a(ii+1)-despl-refay
           az=a(ii+2)-despl-refaz
           radius(i)=sqrt(ax*ax+ay*ay+az*az)
        enddo
     endif

     call indexx(n,radius,indx)  !Indx is sorted by radius measured from
                                 !0,0 when rho<1e12 and from where density
                                 !is higher if rho>1e12.

     RETURN
   END SUBROUTINE
